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


def getPartPlate():
    partPlate = Plate(size=384, well_volume=10, deadspace=0)
    add_volume(partPlate, "A1", 10, bsb.BASIC_SEVA_PARTS["v0.1"]["18"].id)
    add_volume(partPlate, "B1", 10, bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"].id)
    add_volume(partPlate, "C1", 10, bsb.BASIC_SEVA_PARTS["v0.1"]["26"].id)
    partPlate["A1"]["composition"][bsb.BASIC_SEVA_PARTS["v0.1"]
                                   ["18"].id]["concentration"] = 40  # ng / ul
    partPlate["B1"]["composition"][bsb.BASIC_CDS_PARTS["v0.1"]
                                   ["sfGFP"].id]["concentration"] = 30  # ng / ul
    partPlate["C1"]["composition"][bsb.BASIC_SEVA_PARTS["v0.1"]
                                   ["26"].id]["concentration"] = 50  # ng / ul

    return partPlate


def test_echo_instructions_small_build(small_build_example):
    linker_plate = getLinkerPlate()
    part_plate = getPartPlate()
    # print(part_plate["A1"])
    # print(part_plate["B1"])
    # print(part_plate["C1"])
    echo_clips_zippath = bsb.export_echo_clips_instructions(
        small_build_example, linker_plate=linker_plate, part_plate=part_plate)
    assert 1 == 1
