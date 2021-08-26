from .main import BasicBuild
from .csv_export import export_csvs
import zipfile
import os
import pandas as pd
from pathlib import Path
from datetime import datetime
import math
import string
import csv
from functools import reduce
from platemap import Plate, assign_source_wells, find_well, remove_volume
from typing import Literal
from collections import defaultdict


CLIP_VOLUME = 500
BUFFER_VOLUME = 500
TOTAL_VOLUME = 5000


def export_echo_assembly(
    basic_build: BasicBuild,
    path: str = None,
    buffer_well: str = "A1",
    water_well: str = "B1",
    alternate_well: bool = False,
    assemblies_per_clip: int = 20,
    clips_plate_size: Literal[6, 24, 96, 384, 1536] = 384,
    assemblies_plate_size: Literal[6, 24, 96, 384, 1536] = 96,
) -> None:
    """Writes automation scripts for a echo liquid handler to build assemblies from clips.

    Args:
        path (optional): path to zipped folder of csv files. If none defaults to
            working directory with a time stamped name, output csvs is created.
        buffer_well (optional): location in 6 well plate of assembly buffer.
        water_well (optional): location in 6 well plate of dH20.
        alternate_well (optional): specifies whether alternating wells are to be used in the input 384 well plate.
        assemblies_per_clip (optional): number of assemblies each purified clip reaction can support.
        clips_plate_size (optional): specifies the size of the clips plate. Defaults to 384
        assemblies_plate_size (optional): specifiesthe size of the assemblies plates. Defaults to 96


    Returns:
        str: Path of zip file containing echo automation scripts

    Raises:
        ValueError: If water_well or buffer_well is not in ["A1", "B1", "A2", "B2", "A3", "B3"]; if self contains
            96 or more assemblies or if the build requires equal or more than 384 used clip wells for alternate_well(True)
            or 192 for alternate_well(False).
    """

    if water_well not in ["A1", "B1", "A2", "B2", "A3", "B3"]:
        raise ValueError(
            "Water Well location needs to be within the 6 well plate, between A1 - B3"
        )
    if buffer_well not in ["A1", "B1", "A2", "B2", "A3", "B3"]:
        raise ValueError(
            "Assembly Buffer Well location needs to be within the 6 well plate, between A1 - B3"
        )

    calculated_well_volume = assemblies_per_clip * CLIP_VOLUME
    source_plate = Plate(
        size=clips_plate_size, well_volume=calculated_well_volume, deadspace=0
    )
    destination_plate = Plate(size=assemblies_plate_size)

    try:
        assign_source_wells(
            source_plate,
            reduce(
                lambda a, b: {**a, **b},
                list(
                    map(
                        lambda x: {x[0] + 1: len(x[1][1] * CLIP_VOLUME)},
                        enumerate(basic_build.clips_data.items()),
                    )
                ),
            ),
            alternate_wells=alternate_well,
        )

    except:
        raise ValueError(
            """To many clips in the build to be handled by a single 384 
                source plate, considering you alternate_well setting."""
        )

    dd = defaultdict(list)

    for d in list(
        map(
            lambda well_item: {well_item[1][1]["id"]: well_item[1][0]},
            enumerate(
                filter(lambda x: x[1]["total_volume"], source_plate.contents.items())
            ),
        )
    ):
        for key, value in d.items():
            dd[str(key)].append(value)

    clip_sourceplate_mapping = dict(dd)
    assembly_outputplate_mapping = {}

    if path == None:
        now = datetime.now()
        zip_path = (
            Path.cwd() / f"Echo_Instructions_{now.strftime('%d-%m-%Y_%H.%M.%S')}.zip"
        )
    else:
        zip_path = path
    for index, set_of_full_assemblies in enumerate(
        list(
            basic_build.basic_assemblies[x : x + assemblies_plate_size]
            for x in range(0, len(basic_build.basic_assemblies), assemblies_plate_size)
        )
    ):
        with open(
            Path.cwd() / "echo_clips_{}.csv".format(index + 1), "w", newline=""
        ) as f1, open(
            Path.cwd() / "echo_water_buffer_{}.csv".format(index + 1), "w", newline=""
        ) as f2:
            fieldnames = ["Destination Well", "Source Well", "Transfer Volume"]
            thewriter_clips = csv.DictWriter(f1, fieldnames=fieldnames)
            thewriter_clips.writeheader()
            thewriter_water_buffer = csv.DictWriter(f2, fieldnames=fieldnames)
            thewriter_water_buffer.writeheader()
            for assembly_index, assembly in enumerate(set_of_full_assemblies):
                assembly_outputplate_mapping[
                    str((index * assemblies_plate_size) + (assembly_index + 1))
                ] = (str(index) + "-" + str(destination_plate.wells[assembly_index]))
                for clip in [
                    basic_build.unique_clips.index(clip_reaction) + 1
                    for clip_reaction in assembly._clip_reactions
                ]:
                    thewriter_clips.writerow(
                        {
                            "Destination Well": destination_plate.wells[assembly_index],
                            "Source Well": find_well(source_plate, clip, CLIP_VOLUME),
                            "Transfer Volume": CLIP_VOLUME,
                        }
                    )
                    remove_volume(
                        source_plate,
                        find_well(source_plate, clip, CLIP_VOLUME),
                        CLIP_VOLUME,
                    )
                thewriter_water_buffer.writerow(
                    {
                        "Destination Well": destination_plate.wells[assembly_index],
                        "Source Well": buffer_well,
                        "Transfer Volume": BUFFER_VOLUME,
                    }
                )
                thewriter_water_buffer.writerow(
                    {
                        "Destination Well": destination_plate.wells[assembly_index],
                        "Source Well": water_well,
                        "Transfer Volume": TOTAL_VOLUME
                        - BUFFER_VOLUME
                        - CLIP_VOLUME
                        * len(
                            [
                                basic_build.unique_clips.index(clip_reaction)
                                for clip_reaction in assembly._clip_reactions
                            ]
                        ),
                    }
                )
    csv_zip = export_csvs(
        basic_build, None, clip_sourceplate_mapping, assembly_outputplate_mapping
    )
    with zipfile.ZipFile(csv_zip, "r") as zip_ref:
        zip_ref.extractall()
    with zipfile.ZipFile(zip_path, "w") as my_zip:
        my_zip.write("clips.csv")
        my_zip.write("assemblies.csv")
        os.remove("clips.csv")
        os.remove("assemblies.csv")
        os.remove(csv_zip)
        for file in os.listdir(Path.cwd()):
            if file.startswith("echo_") and file.endswith(".csv"):
                my_zip.write(file)
                os.remove(file)
    return zip_path
