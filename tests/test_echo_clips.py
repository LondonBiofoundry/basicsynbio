from platemap.PlateUtils import add_volume
from platemap.plate import Plate
import basicsynbio as bsb
import pytest
from .test_fixtures import small_build_example


def getLinkerPlate():
    linkerPlate = Plate(size=384, well_volume=10, deadspace=0)
    add_volume(linkerPlate, "A1", 10, "LMP-P")
    add_volume(linkerPlate, "B1", 10, "LMP-S")
    add_volume(linkerPlate, "C1", 10, "LMS-P")
    add_volume(linkerPlate, "D1", 10, "LMS-S")
    return linkerPlate


def test_echo_instructions_small_build(small_build_example):
    linker_plate = getLinkerPlate()
    print(linker_plate["A1"])
    print(linker_plate["B1"])
    print(linker_plate["C1"])
    print(linker_plate["D1"])
    echo_clips_zippath = bsb.export_echo_clips_instructions(
        small_build_example, linker_plate=linker_plate)
    assert 1 == 1
