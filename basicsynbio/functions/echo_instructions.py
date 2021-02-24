from basicsynbio.cam import BasicBuild
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

CLIP_VOLUME = 500
BUFFER_VOLUME = 500
TOTAL_VOLUME = 5000


def export_echo_assembly(
    basic_build: BasicBuild,
    path: str = None,
    buffer_well: str = "A1",
    water_well: str = "B1",
    alternate_well: bool = False,
) -> None:
    """Writes automation scripts for a echo liquid handler to build assemblies from clips.

    Args:
        path (optional): path to zipped folder of csv files. If none defaults to
            working directory with a time stamped name, output csvs is created.
        buffer_well (optional): location in 6 well plate of assembly buffer.
        water_well (optional): location in 6 well plate of dH20.
        alternate_well (optional): specifies whether alternating wells are to be used in the input 384 well plate.

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

    source_plate = Plate(size=384, well_volume=40000, deadspace=20000)
    destination_plate = Plate(size=96)
    try:
        assign_source_wells(
            source_plate,
            reduce(
                lambda a, b: {**a, **b},
                list(
                    map(
                        lambda x: {x[0]: len(x[1][1] * CLIP_VOLUME)},
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

    if path == None:
        now = datetime.now()
        zip_path = (
            Path.cwd() / f"Echo_Instructions_{now.strftime('%d-%m-%Y_%H.%M.%S')}.zip"
        )
    else:
        zip_path = path
    for index, set_of_96_assemblies in enumerate(
        list(
            basic_build.basic_assemblies[x : x + 96]
            for x in range(0, len(basic_build.basic_assemblies), 96)
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
            for index, assembly in enumerate(set_of_96_assemblies):
                for clip in [
                    basic_build.unique_clips.index(clip_reaction)
                    for clip_reaction in assembly.clip_reactions
                ]:
                    thewriter_clips.writerow(
                        {
                            "Destination Well": destination_plate.wells[index],
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
                        "Destination Well": destination_plate.wells[index],
                        "Source Well": buffer_well,
                        "Transfer Volume": BUFFER_VOLUME,
                    }
                )
                thewriter_water_buffer.writerow(
                    {
                        "Destination Well": destination_plate.wells[index],
                        "Source Well": water_well,
                        "Transfer Volume": TOTAL_VOLUME
                        - BUFFER_VOLUME
                        - CLIP_VOLUME
                        * len(
                            [
                                basic_build.unique_clips.index(clip_reaction)
                                for clip_reaction in assembly.clip_reactions
                            ]
                        ),
                    }
                )
    with zipfile.ZipFile(zip_path, "w") as my_zip:
        for file in os.listdir(Path.cwd()):
            if file.startswith("echo_") and file.endswith(".csv"):
                my_zip.write(file)
                os.remove(file)
    return zip_path
