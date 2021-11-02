import csv
from datetime import datetime
from .main import BasicBuild
import os
import pandas as pd
from pathlib import Path
from typing import Dict, List
import zipfile


def export_csvs(
    basic_build: BasicBuild,
    path: str = None,
    clip_plate_mapping: Dict[str, List[str]] = None,
    assembly_plate_mapping: Dict[str, str] = None,
):
    """Writes information about each clip_data and assembly to
    two dependent CSV files in the same folder the command
    is executed.

    Args:
        path (optional): path to zipped folder of csv files. If none defaults to
            working directory with a time stamped name, output csvs is created.
        clip_plate_mapping (optional): A dictionary with keys containing clip
            indicies and values containing lists of well locations's where the
            clips are stored.
        assembly_plate_mapping (optional): A dictionary with keys containing
            assembly indicies and values containing the assemblies location plate
            followed by a dash then the well locations for example {1:'0-A1',...}


    Returns:
        str: filepath of created zip containing the CSV files
    """
    if path == None:
        now = datetime.now()
        zip_path = Path.cwd() / f"build_{now.strftime('%d-%m-%Y_%H.%M.%S')}.zip"
    else:
        zip_path = path
    with open(Path.cwd() / "clips.csv", "w", newline="") as f:
        fieldnames = [
            "Clip Index",
            "Prefix ID",
            "Part ID",
            "Part Name",
            "Part mass for 30 μL clip reaction (ng)",
            "Suffix ID",
            "Total assemblies",
            "Assembly indexes",
            "Clip plate mapping",
        ]
        thewriter = csv.DictWriter(f, fieldnames=fieldnames)
        thewriter.writeheader()
        for index, clip_data in enumerate(basic_build.clips_data.items()):
            thewriter.writerow(
                {
                    "Clip Index": index + 1,
                    "Prefix ID": clip_data[0]._prefix.prefix_id,
                    "Part ID": clip_data[0]._part.id,
                    "Part Name": clip_data[0]._part.name,
                    "Part mass for 30 μL clip reaction (ng)": clip_data[
                        0
                    ]._part.clip_mass(),
                    "Suffix ID": clip_data[0]._suffix.suffix_id,
                    "Total assemblies": len(clip_data[1]),
                    "Assembly indexes": [
                        basic_build.basic_assemblies.index(assembly) + 1
                        for assembly in clip_data[1]
                    ],
                    "Clip plate mapping": clip_plate_mapping[str(index + 1)]
                    if clip_plate_mapping
                    else "N/A",
                }
            )
    with open(Path.cwd() / "assemblies.csv", "w", newline="") as f:
        fieldnames = [
            "Assembly Index",
            "Assembly ID",
            "Clip indexes",
            "Assembly plate mapping",
        ]
        thewriter = csv.DictWriter(f, fieldnames=fieldnames)
        thewriter.writeheader()
        for index, assembly in enumerate(basic_build.basic_assemblies):
            thewriter.writerow(
                {
                    "Assembly Index": index + 1,
                    "Assembly ID": assembly.id,
                    "Clip indexes": [
                        basic_build.unique_clips.index(clip_reaction) + 1
                        for clip_reaction in assembly._clip_reactions
                    ],
                    "Assembly plate mapping": assembly_plate_mapping[str(index + 1)]
                    if assembly_plate_mapping
                    else "N/A",
                }
            )
    with zipfile.ZipFile(zip_path, "w") as my_zip:
        try:
            my_zip.write("assemblies.csv")
            my_zip.write("clips.csv")
        finally:
            my_zip.close()
    os.remove(Path.cwd() / "assemblies.csv")
    os.remove(Path.cwd() / "clips.csv")
    return zip_path
