from platemap.PlateUtils import add_volume
from platemap.plate import Plate
import basicsynbio as bsb
import zipfile
import os
import pandas as pd
import numpy as np
from pathlib import Path
import pytest
from .test_fixtures import small_build_example


def getLinkerPlate():
    linkerPlate = Plate(size=384, well_volume=10, deadspace=0)
    add_volume(linkerPlate, "A1", 10, "LMP-P")
    add_volume(linkerPlate, "B1", 10, "LMP-S")
    add_volume(linkerPlate, "C1", 10, "LMS-P")
    add_volume(linkerPlate, "D1", 10, "LMS-S")
    return linkerPlate


def getPartPlate():
    partPlate = Plate(size=384, well_volume=10, deadspace=0)
    add_volume(partPlate, "A1", 10, bsb.BASIC_SEVA_PARTS["v0.1"]["18"].id)
    add_volume(partPlate, "B1", 10, bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"].id)
    add_volume(partPlate, "C1", 10, bsb.BASIC_SEVA_PARTS["v0.1"]["26"].id)
    partPlate["A1"]["composition"][bsb.BASIC_SEVA_PARTS["v0.1"]["18"].id][
        "concentration"
    ] = 40  # ng / ul
    partPlate["B1"]["composition"][bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"].id][
        "concentration"
    ] = 30  # ng / ul
    partPlate["C1"]["composition"][bsb.BASIC_SEVA_PARTS["v0.1"]["26"].id][
        "concentration"
    ] = 50  # ng / ul

    return partPlate


def test_echo_instructions_small_build(small_build_example):
    linker_plate = getLinkerPlate()
    part_plate = getPartPlate()
    echo_clips_zippath = bsb.export_echo_clips_instructions(
        small_build_example, linker_plate=linker_plate, part_plate=part_plate
    )
    with zipfile.ZipFile(echo_clips_zippath, "r") as zip_ref:
        try:
            zip_ref.extractall()
        finally:
            zip_ref.close()
            os.remove(echo_clips_zippath)
    stage1 = pd.read_csv(Path.cwd() / "stage_1_half_linkers.csv")
    stage2 = pd.read_csv(Path.cwd() / "stage_2_parts.csv")
    stage3 = pd.read_csv(Path.cwd() / "stage_3_water_buffer.csv")
    os.remove(Path.cwd() / "stage_1_half_linkers.csv")
    os.remove(Path.cwd() / "stage_2_parts.csv")
    os.remove(Path.cwd() / "stage_3_water_buffer.csv")
    expected_stage1 = [
        ["A1", "C1", 0.7],
        ["A1", "B1", 0.7],
        ["B1", "A1", 0.7],
        ["B1", "D1", 0.7],
        ["C1", "C1", 0.7],
        ["C1", "B1", 0.7],
    ]
    expected_stage2 = [
        ["A1", "A1", 2.7],
        ["B1", "B1", 3.4],
        ["C1", "C1", 2.0],
    ]
    expected_stage3 = [
        ["A1", "A1", 6.7],
        ["A1", "B1", 9.2],
        ["B1", "A1", 6.7],
        ["B1", "B1", 8.5],
        ["C1", "A1", 6.7],
        ["C1", "B1", 9.9],
    ]
    assert expected_stage1 == stage1.to_numpy().tolist()
    assert expected_stage2 == stage2.to_numpy().tolist()
    assert expected_stage3 == stage3.to_numpy().tolist()


def test_echo_instructions_small_build_default_plate(small_build_example):
    part_plate = getPartPlate()
    echo_clips_zippath = bsb.export_echo_clips_instructions(
        small_build_example, part_plate=part_plate
    )
    with zipfile.ZipFile(echo_clips_zippath, "r") as zip_ref:
        try:
            zip_ref.extractall()
        finally:
            zip_ref.close()
            os.remove(echo_clips_zippath)
    stage1 = pd.read_csv(Path.cwd() / "stage_1_half_linkers.csv")
    stage2 = pd.read_csv(Path.cwd() / "stage_2_parts.csv")
    stage3 = pd.read_csv(Path.cwd() / "stage_3_water_buffer.csv")
    os.remove(Path.cwd() / "stage_1_half_linkers.csv")
    os.remove(Path.cwd() / "stage_2_parts.csv")
    os.remove(Path.cwd() / "stage_3_water_buffer.csv")
    expected_stage1 = [
        ["A1", "C15", 0.7],
        ["A1", "A13", 0.7],
        ["B1", "C13", 0.7],
        ["B1", "A15", 0.7],
        ["C1", "C15", 0.7],
        ["C1", "A13", 0.7],
    ]
    expected_stage2 = [
        ["A1", "A1", 2.7],
        ["B1", "B1", 3.4],
        ["C1", "C1", 2.0],
    ]
    expected_stage3 = [
        ["A1", "A1", 6.7],
        ["A1", "B1", 9.2],
        ["B1", "A1", 6.7],
        ["B1", "B1", 8.5],
        ["C1", "A1", 6.7],
        ["C1", "B1", 9.9],
    ]
    assert expected_stage1 == stage1.to_numpy().tolist()
    assert expected_stage2 == stage2.to_numpy().tolist()
    assert expected_stage3 == stage3.to_numpy().tolist()
