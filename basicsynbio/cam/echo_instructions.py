from .default_linker_plate import LINKER_384_PLATE
from .main import BasicBuild, BuildEncoder
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
import json


CLIP_VOLUME = 500
BUFFER_VOLUME = 500
TOTAL_VOLUME = 5000


def export_echo_assembly_instructions(
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
    my_zip.close()
    return zip_path


def export_echo_clips_instructions(
    basic_build: BasicBuild,
    part_plate: Plate,
    linker_plate: Plate = LINKER_384_PLATE,
    fold_dilution: float = 0.7,
    buffer_well: str = "A1",
    water_well: str = "B1",
    path: str = None,
) -> None:
    """Export automation instuctions to build the basic clips present within the basic build.

    A function to export the necessary csvs files for the Labcyte Echo liquid handling robot to transform half
    linkers and parts into basic clips ready for assembly into full basic assemblies.

    The assembly occurs in three stages each with their corresponding csv file, the source plates for each stage
    change however the desination plates do not. The three stages are:

        1. Adding the half-linkers to the desination plate
        2. Adding the parts to the desination plate
        3. Adding water and the buffer to the desination plate

    The input parameter, linker_plate, is the plate containing the half-linkers, the source plate supplied for the
    first stage. The input parameter, parts_plate, is the plate containing the parts, the source plate supplied for
    the second stage. The input parameters buffer_well and water_well define where in the source plate for the first
    stage the buffer and water are located.

    Args:
        basic_build: the basic build to export the instructions for.
        parts_plate: A platemap instance containing the half linker locations and volumes.
        linker_plate: A platemap instance containing the half linker locations and volumes.
        fold_dilution: The manual workflow to sythesise basic clips from half linkers and parts forms 30 µl of
            each linker. The Labcyte Echo robot is designed to have the desination plate upside down and volumes
            over 20µl have a risk of being lost. As the lab instructions generated in this function target this
            machine we default to a fold dilution of 0.7. to dilute 30µl to approximately 20µl. We obtain exactly
            20µl by adjusting the amount of water we add as the final reagent added to the destination plate.
        buffer_well: The location of the buffer well on the source plate.
        water_well: The location of the water well on the source plate.
        path (optional): path to zipped folder of csv files. If none defaults to working directory with a time
            stamped name, output csvs is created.


    Returns:
        str: The Path of zip file containing the Labcyte Echo clips automation scripts

    Raises:
        ValueError: Export clips is currently limited to 96 clips, please reduce the number of clips.
    """
    print("Planning destination plate")
    if len(basic_build.unique_clips) > 96:
        raise ValueError(
            "Export clips is currently limited to 96 clips, please reduce the number of clips"
        )
    if water_well not in ["A1", "B1", "A2", "B2", "A3", "B3"]:
        raise ValueError(
            "Water Well location needs to be within the 6 well plate, between A1 - B3"
        )
    if buffer_well not in ["A1", "B1", "A2", "B2", "A3", "B3"]:
        raise ValueError(
            "Assembly Buffer Well location needs to be within the 6 well plate, between A1 - B3"
        )
    # Defining function variables
    HALF_LINKER_VOLUME = 1 * fold_dilution
    destination_plate = Plate(size=96, well_volume=20)

    # Stage 1
    stage_1_liquid_transfers = []
    stage_2_liquid_transfers = []
    stage_3_liquid_transfers = []

    for index, clip in enumerate(basic_build.unique_clips):
        # Defining the location of each clip in the desination plate.
        destination_plate.set_well_id(
            destination_plate.wells[index], "CR{}".format(index + 1)
        )
        # Finding the location of the half-linker in the linker plate
        prefix_half_linker_id, suffix_half_linker_id = clip.linker_half_ids()
        prefix_well = find_well(linker_plate, prefix_half_linker_id, HALF_LINKER_VOLUME)
        if prefix_well == 0:  # The value returned in no well was found
            raise ValueError(
                "The half linker {} is not in the source plate".format(
                    prefix_half_linker_id
                )
            )
        suffix_well = find_well(linker_plate, suffix_half_linker_id, HALF_LINKER_VOLUME)
        if suffix_well == 0:  # The value returned in no well was found
            raise ValueError(
                "The half linker {} is not in the source plate".format(
                    suffix_half_linker_id
                )
            )
        # Adding the transfer instructions to the list
        # Adding prefix half linkers
        stage_1_liquid_transfers.append(
            {
                "Destination Well": destination_plate.wells[index],
                "Source Well": prefix_well,
                "Transfer Volume": HALF_LINKER_VOLUME,
            }
        )
        # Adding suffix half linkers
        stage_1_liquid_transfers.append(
            {
                "Destination Well": destination_plate.wells[index],
                "Source Well": suffix_well,
                "Transfer Volume": HALF_LINKER_VOLUME,
            }
        )

        # Stage 2
        basic_part = list(clip.clip_items())[1]
        # Calculate volume required of part
        required_mass_nano_grams = basic_part.clip_mass(clip_vol=30 * fold_dilution)
        part_well = find_well(part_plate, basic_part.id, 0)
        if part_well == 0:  # The value returned in no well was found
            raise ValueError(
                "The part {} is not in the part plate".format(basic_part.name)
            )
        try:
            required_volume = (
                required_mass_nano_grams
                / part_plate[part_well]["composition"][basic_part.id]["concentration"]
            )
        except:
            raise ValueError(
                "The part {} in the part_plate is does not have a concentration defined, please define and retry".format(
                    basic_part.name
                )
            )
        required_volume_1dp = round(required_volume, 1)

        # Finding the location of the parts in the part plate
        part_well_with_requiured_volume = find_well(
            part_plate, basic_part.id, required_volume_1dp
        )
        if (
            part_well_with_requiured_volume == 0
        ):  # The value returned in no well was found
            raise ValueError(
                "The part {} is not in the part plate".format(basic_part.name)
            )

        # Adding the transfer instructions to the list
        stage_2_liquid_transfers.append(
            {
                "Destination Well": destination_plate.wells[index],
                "Source Well": part_well_with_requiured_volume,
                "Transfer Volume": required_volume_1dp,
            }
        )

        # Stage 3
        # Add buffer
        stage_3_liquid_transfers.append(
            {
                "Destination Well": destination_plate.wells[index],
                "Source Well": buffer_well,
                "Transfer Volume": round(20 / 3, 1),
            }
        )
        # Add Water
        water_volume = round(
            20 - round(20 / 3, 1) - required_volume_1dp - (2 * HALF_LINKER_VOLUME), 1
        )
        if water_volume < 0:
            raise ValueError(
                "Cannot add more water than the destination plate can hold, increase the concentration of the parts"
            )
        stage_3_liquid_transfers.append(
            {
                "Destination Well": destination_plate.wells[index],
                "Source Well": water_well,
                "Transfer Volume": water_volume,
            }
        )

    # Write transfer steps to CSV
    if path == None:
        now = datetime.now()
        zip_path = (
            Path.cwd() / f"Echo_Instructions_{now.strftime('%d-%m-%Y_%H.%M.%S')}.zip"
        )
    else:
        zip_path = path

    with open(Path.cwd() / "stage_1_half_linkers.csv", "w", newline="") as f1, open(
        Path.cwd() / "stage_2_parts.csv", "w", newline=""
    ) as f2, open(Path.cwd() / "stage_3_water_buffer.csv", "w", newline="") as f3:
        fieldnames = ["Destination Well", "Source Well", "Transfer Volume"]
        w1 = csv.DictWriter(f1, fieldnames=fieldnames)
        w1.writeheader()
        w2 = csv.DictWriter(f2, fieldnames=fieldnames)
        w2.writeheader()
        w3 = csv.DictWriter(f3, fieldnames=fieldnames)
        w3.writeheader()
        for transfer in stage_1_liquid_transfers:
            w1.writerow(transfer)
        for transfer in stage_2_liquid_transfers:
            w2.writerow(transfer)
        for transfer in stage_3_liquid_transfers:
            w3.writerow(transfer)

    with zipfile.ZipFile(zip_path, "w") as my_zip:
        my_zip.write("stage_1_half_linkers.csv")
        my_zip.write("stage_2_parts.csv")
        my_zip.write("stage_3_water_buffer.csv")
        os.remove("stage_1_half_linkers.csv")
        os.remove("stage_2_parts.csv")
        os.remove("stage_3_water_buffer.csv")

    my_zip.close()
    # print("Stage 3")
    for transfer in stage_3_liquid_transfers:
        print(transfer)
    return zip_path
