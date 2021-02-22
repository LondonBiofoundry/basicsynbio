import basicsynbio as bsb
import pytest

from test_fixtures import (
    gfp_basicpart,
    gfp_seqrec,
    gfp_orf_seq,
    cmr_p15a_basicpart,
    cmr_p15a_backbone,
    five_part_assembly_parts,
    five_part_assembly_linkers,
    five_part_assembly,
    gfp_orf_seqrec,
    gfp_orf_basicpart,
    bseva_68_seqrec,
    bsai_part_seqrec,
    promoter_assemblies_build,
    promoter_assemblies_json,
    gfp_part_final_conc,
    small_build_example,
    all_promoter_assemblies_build,
)


def test_echo_instructions_small_build(small_build_example):
    import zipfile
    import os
    import pandas as pd
    import numpy as np
    from pathlib import Path

    echozippath = bsb.export_echo_assembly(small_build_example)
    with zipfile.ZipFile(echozippath, "r") as zip_ref:
        try:
            zip_ref.extractall("ECHO_CSVS")
        finally:
            zip_ref.close()
            os.remove(echozippath)
    echo_clips = pd.read_csv(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    echo_water_buffer = pd.read_csv(
        Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv"
    )
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv")
    os.rmdir(Path.cwd() / "ECHO_CSVS")
    expected_clips = [
        ["A1", "A1", 500],
        ["A1", "B1", 500],
        ["B1", "C1", 500],
        ["B1", "B1", 500],
    ]
    expected_water_buffer = [
        ["A1", "A1", 500],
        ["A1", "B1", 3500],
        ["B1", "A1", 500],
        ["B1", "B1", 3500],
    ]
    assert expected_clips == echo_clips.to_numpy().tolist()
    assert expected_water_buffer == echo_water_buffer.to_numpy().tolist()


def test_echo_instructions_small_build_useAllWell_False(small_build_example):
    import zipfile
    import os
    import pandas as pd
    import numpy as np
    from pathlib import Path

    echozippath = bsb.export_echo_assembly(small_build_example, alternate_well=True)
    with zipfile.ZipFile(echozippath, "r") as zip_ref:
        try:
            zip_ref.extractall("ECHO_CSVS")
        finally:
            zip_ref.close()
            os.remove(echozippath)
    echo_clips = pd.read_csv(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    echo_water_buffer = pd.read_csv(
        Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv"
    )
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv")
    os.rmdir(Path.cwd() / "ECHO_CSVS")
    expected_clips = [
        ["A1", "A1", 500],
        ["A1", "C1", 500],
        ["B1", "E1", 500],
        ["B1", "C1", 500],
    ]
    expected_water_buffer = [
        ["A1", "A1", 500],
        ["A1", "B1", 3500],
        ["B1", "A1", 500],
        ["B1", "B1", 3500],
    ]
    assert expected_clips == echo_clips.to_numpy().tolist()
    assert expected_water_buffer == echo_water_buffer.to_numpy().tolist()


def test_echo_instruction_assert_buffer_water_well_errors(small_build_example):
    with pytest.raises(
        ValueError,
        match="Water Well location needs to be within the 6 well plate, between A1 - B3",
    ):
        echozippath = bsb.export_echo_assembly(small_build_example, waterWell="D1")
    with pytest.raises(
        ValueError,
        match="Assembly Buffer Well location needs to be within the 6 well plate, between A1 - B3",
    ):
        echozippath = bsb.export_echo_assembly(small_build_example, bufferWell="D1")


def test_echo_instruction_assert_buffer_water_well_errors(promoter_assemblies_build):
    # This Build contains ~180 assemblies function should raise errors for builds
    # with more than 96 assemblies
    with pytest.raises(ValueError, match=r".*To many assemblies in the build.*"):
        echozippath = bsb.export_echo_assembly(promoter_assemblies_build)


def test_echo_path(small_build_example):
    from pathlib import Path
    import os

    zippath = Path.cwd() / "uniquezippath.zip"
    realzippath = bsb.export_echo_assembly(small_build_example, path=zippath)
    os.remove(zippath)
    assert zippath == realzippath
