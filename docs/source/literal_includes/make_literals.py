# -*- coding: utf-8 -*-
"""Make files to include in sphinx docs.

Running this module generates additional files in this directory, 
which are then included in the documentation. This is useful for 
ensuring documented variables are current.

"""

import basicsynbio as bsb
import zipfile
from pathlib import Path
import os

def build_json():
    promoter_assemblies = (bsb.BasicAssembly(
        f"promoter_construct_{ind}",
        bsb.BASIC_SEVA_PARTS["v0.1"]["26"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
        promoter,
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR1-RBS2"],
        bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"]
        ) for ind, promoter in enumerate(bsb.BASIC_PROMOTER_PARTS["v0.1"].values()))
    build = bsb.BasicBuild(*promoter_assemblies)
    return build

def export_csvs(build):
    build.export_csvs("build_csvs.zip")

def export_json(build):
    import json
    with open("build.json", "w") as json_file:
        json.dump(build, json_file, cls=bsb.BuildEncoder, indent=4, ensure_ascii=False)

if __name__ == "__main__":
    export_json(build_json())
    BUILD_ZIP_PATH = "build_csvs.zip"
    try:
        export_csvs(build_json())
        with zipfile.ZipFile(BUILD_ZIP_PATH, 'r') as myzip:
            myzip.extractall()
    finally:
        os.remove(BUILD_ZIP_PATH)


    