import basicsynbio as bsb
import pytest

from .test_fixtures import (
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
    promoter_assemblies_build_more_than_384,
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
        Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_1.csv"
    )
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_1.csv")
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
        Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_1.csv"
    )
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_1.csv")
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
        echozippath = bsb.export_echo_assembly(small_build_example, water_well="D1")
    with pytest.raises(
        ValueError,
        match="Assembly Buffer Well location needs to be within the 6 well plate, between A1 - B3",
    ):
        echozippath = bsb.export_echo_assembly(small_build_example, buffer_well="D1")


def test_echo_path(small_build_example):
    from pathlib import Path
    import os

    zippath = Path.cwd() / "uniquezippath.zip"
    realzippath = bsb.export_echo_assembly(small_build_example, path=zippath)
    os.remove(zippath)
    assert zippath == realzippath


def test_multiple_files_made_more_than_96_assemblies(promoter_assemblies_build):
    import zipfile
    import os
    import pandas as pd
    import numpy as np
    from pathlib import Path

    echozippath = bsb.export_echo_assembly(promoter_assemblies_build)
    with zipfile.ZipFile(echozippath, "r") as zip_ref:
        try:
            zip_ref.extractall("ECHO_CSVS")
        finally:
            zip_ref.close()
            os.remove(echozippath)
    echo_clips = pd.read_csv(Path.cwd() / "ECHO_CSVS" / "echo_clips_2.csv")
    echo_water_buffer_2 = pd.read_csv(
        Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv"
    )
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_2.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv")
    os.rmdir(Path.cwd() / "ECHO_CSVS")
    expected_water_buffer_2 = [
        ["A1", "A1", 500],
        ["A1", "B1", 3000],
        ["B1", "A1", 500],
        ["B1", "B1", 3000],
        ["C1", "A1", 500],
        ["C1", "B1", 3000],
        ["D1", "A1", 500],
        ["D1", "B1", 3000],
        ["E1", "A1", 500],
        ["E1", "B1", 3000],
        ["F1", "A1", 500],
        ["F1", "B1", 3000],
        ["G1", "A1", 500],
        ["G1", "B1", 3000],
        ["H1", "A1", 500],
        ["H1", "B1", 3000],
        ["A2", "A1", 500],
        ["A2", "B1", 3000],
        ["B2", "A1", 500],
        ["B2", "B1", 3000],
        ["C2", "A1", 500],
        ["C2", "B1", 3000],
        ["D2", "A1", 500],
        ["D2", "B1", 3000],
        ["E2", "A1", 500],
        ["E2", "B1", 3000],
        ["F2", "A1", 500],
        ["F2", "B1", 3000],
        ["G2", "A1", 500],
        ["G2", "B1", 3000],
        ["H2", "A1", 500],
        ["H2", "B1", 3000],
        ["A3", "A1", 500],
        ["A3", "B1", 3000],
        ["B3", "A1", 500],
        ["B3", "B1", 3000],
        ["C3", "A1", 500],
        ["C3", "B1", 3000],
        ["D3", "A1", 500],
        ["D3", "B1", 3000],
        ["E3", "A1", 500],
        ["E3", "B1", 3000],
        ["F3", "A1", 500],
        ["F3", "B1", 3000],
        ["G3", "A1", 500],
        ["G3", "B1", 3000],
        ["H3", "A1", 500],
        ["H3", "B1", 3000],
        ["A4", "A1", 500],
        ["A4", "B1", 3000],
        ["B4", "A1", 500],
        ["B4", "B1", 3000],
        ["C4", "A1", 500],
        ["C4", "B1", 3000],
        ["D4", "A1", 500],
        ["D4", "B1", 3000],
        ["E4", "A1", 500],
        ["E4", "B1", 3000],
        ["F4", "A1", 500],
        ["F4", "B1", 3000],
        ["G4", "A1", 500],
        ["G4", "B1", 3000],
        ["H4", "A1", 500],
        ["H4", "B1", 3000],
        ["A5", "A1", 500],
        ["A5", "B1", 3000],
        ["B5", "A1", 500],
        ["B5", "B1", 3000],
        ["C5", "A1", 500],
        ["C5", "B1", 3000],
        ["D5", "A1", 500],
        ["D5", "B1", 3000],
        ["E5", "A1", 500],
        ["E5", "B1", 3000],
        ["F5", "A1", 500],
        ["F5", "B1", 3000],
        ["G5", "A1", 500],
        ["G5", "B1", 3000],
        ["H5", "A1", 500],
        ["H5", "B1", 3000],
        ["A6", "A1", 500],
        ["A6", "B1", 3000],
        ["B6", "A1", 500],
        ["B6", "B1", 3000],
        ["C6", "A1", 500],
        ["C6", "B1", 3000],
        ["D6", "A1", 500],
        ["D6", "B1", 3000],
        ["E6", "A1", 500],
        ["E6", "B1", 3000],
        ["F6", "A1", 500],
        ["F6", "B1", 3000],
        ["G6", "A1", 500],
        ["G6", "B1", 3000],
        ["H6", "A1", 500],
        ["H6", "B1", 3000],
        ["A7", "A1", 500],
        ["A7", "B1", 3000],
        ["B7", "A1", 500],
        ["B7", "B1", 3000],
        ["C7", "A1", 500],
        ["C7", "B1", 3000],
        ["D7", "A1", 500],
        ["D7", "B1", 3000],
        ["E7", "A1", 500],
        ["E7", "B1", 3000],
        ["F7", "A1", 500],
        ["F7", "B1", 3000],
        ["G7", "A1", 500],
        ["G7", "B1", 3000],
        ["H7", "A1", 500],
        ["H7", "B1", 3000],
        ["A8", "A1", 500],
        ["A8", "B1", 3000],
        ["B8", "A1", 500],
        ["B8", "B1", 3000],
        ["C8", "A1", 500],
        ["C8", "B1", 3000],
        ["D8", "A1", 500],
        ["D8", "B1", 3000],
        ["E8", "A1", 500],
        ["E8", "B1", 3000],
        ["F8", "A1", 500],
        ["F8", "B1", 3000],
        ["G8", "A1", 500],
        ["G8", "B1", 3000],
        ["H8", "A1", 500],
        ["H8", "B1", 3000],
        ["A9", "A1", 500],
        ["A9", "B1", 3000],
        ["B9", "A1", 500],
        ["B9", "B1", 3000],
        ["C9", "A1", 500],
        ["C9", "B1", 3000],
        ["D9", "A1", 500],
        ["D9", "B1", 3000],
        ["E9", "A1", 500],
        ["E9", "B1", 3000],
        ["F9", "A1", 500],
        ["F9", "B1", 3000],
        ["G9", "A1", 500],
        ["G9", "B1", 3000],
        ["H9", "A1", 500],
        ["H9", "B1", 3000],
        ["A10", "A1", 500],
        ["A10", "B1", 3000],
        ["B10", "A1", 500],
        ["B10", "B1", 3000],
        ["C10", "A1", 500],
        ["C10", "B1", 3000],
        ["D10", "A1", 500],
        ["D10", "B1", 3000],
        ["E10", "A1", 500],
        ["E10", "B1", 3000],
        ["F10", "A1", 500],
        ["F10", "B1", 3000],
        ["G10", "A1", 500],
        ["G10", "B1", 3000],
        ["H10", "A1", 500],
        ["H10", "B1", 3000],
        ["A11", "A1", 500],
        ["A11", "B1", 3000],
        ["B11", "A1", 500],
        ["B11", "B1", 3000],
        ["C11", "A1", 500],
        ["C11", "B1", 3000],
        ["D11", "A1", 500],
        ["D11", "B1", 3000],
        ["E11", "A1", 500],
        ["E11", "B1", 3000],
        ["F11", "A1", 500],
        ["F11", "B1", 3000],
        ["G11", "A1", 500],
        ["G11", "B1", 3000],
        ["H11", "A1", 500],
        ["H11", "B1", 3000],
        ["A12", "A1", 500],
        ["A12", "B1", 3000],
        ["B12", "A1", 500],
        ["B12", "B1", 3000],
        ["C12", "A1", 500],
        ["C12", "B1", 3000],
        ["D12", "A1", 500],
        ["D12", "B1", 3000],
        ["E12", "A1", 500],
        ["E12", "B1", 3000],
    ]
    assert expected_water_buffer_2 == echo_water_buffer_2.to_numpy().tolist()


def test_echo_instructions_too_many_clips(promoter_assemblies_build_more_than_384):
    import zipfile
    import os
    import pandas as pd
    import numpy as np
    from pathlib import Path

    with pytest.raises(ValueError):
        bsb.export_echo_assembly(promoter_assemblies_build_more_than_384)
