from .main import BasicBuild
import zipfile
import os
import pandas as pd
from pathlib import Path
from datetime import datetime
import csv
    
def export_csvs(basic_build: BasicBuild, path: str = None):
    """Writes information about each clip_data and assembly to
    two dependent CSV files in the same folder the command
    is executed.

    Args:
        path (optional): path to zipped folder of csv files. If none defaults to
            working directory with a time stamped name, output csvs is created.

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
            "Part suggested stock concentration (ng/µL)",
            "Part stock per 30 µL clip (µL)",
            "Suffix ID",
            "Total assemblies",
            "Assembly indexes",
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
                    "Part suggested stock concentration (ng/µL)": clip_data[
                        0
                    ]._part.concentration(),
                    "Part stock per 30 µL clip (µL)": 1,
                    "Suffix ID": clip_data[0]._suffix.suffix_id,
                    "Total assemblies": len(clip_data[1]),
                    "Assembly indexes": [
                        basic_build.basic_assemblies.index(assembly) + 1
                        for assembly in clip_data[1]
                    ],
                }
            )
    with open(Path.cwd() / "assemblies.csv", "w", newline="") as f:
        fieldnames = ["Assembly Index", "Assembly ID", "Clip indexes"]
        thewriter = csv.DictWriter(f, fieldnames=fieldnames)
        thewriter.writeheader()
        for index, assembly in enumerate(basic_build.basic_assemblies):
            thewriter.writerow(
                {
                    "Assembly Index": index + 1,
                    "Assembly ID": assembly.id,
                    "Clip indexes": [
                        basic_build.unique_clips.index(clip_reaction) + 1
                        for clip_reaction in assembly.clip_reactions
                    ],
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