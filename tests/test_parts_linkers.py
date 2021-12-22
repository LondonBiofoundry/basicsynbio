import basicsynbio as bsb
from .test_fixtures import bseva_68_seqrec
from .functions import compare_seqrec_instances


def test_bseva_dict(bseva_68_seqrec):
    print(bsb.BASIC_SEVA_PARTS["v0.1"].keys())
    bseva_68_part = bsb.BASIC_SEVA_PARTS["v0.1"]["68"]
    assert compare_seqrec_instances(bseva_68_part, bseva_68_seqrec) == True


def test_bpromoter_dict():
    from Bio import SeqIO

    bpromoters_handle = "./basicsynbio/parts_linkers/BASIC_promoter_collection_v01.gb"
    bpromoter_seqrecs = SeqIO.parse(bpromoters_handle, "genbank")
    for seqrec in bpromoter_seqrecs:
        collection_key = seqrec.id
        setattr(seqrec, "id", bsb.cam.seqrecord_hexdigest(seqrec))
        assert (
            compare_seqrec_instances(
                bsb.BASIC_PROMOTER_PARTS["v0.1"][collection_key], seqrec
            )
            == True
        )


def test_bcds_dict():
    from Bio import SeqIO

    bcds_handle = "./basicsynbio/parts_linkers/BASIC_CDS_collection_v01.gb"
    bcds_seqrecs = SeqIO.parse(bcds_handle, "genbank")
    for seqrec in bcds_seqrecs:
        collection_key = seqrec.id
        setattr(seqrec, "id", bsb.cam.seqrecord_hexdigest(seqrec))
        assert (
            compare_seqrec_instances(
                bsb.BASIC_CDS_PARTS["v0.1"][collection_key], seqrec
            )
            == True
        )


def test_basic_cds_parts_v02():
    assert len(bsb.BASIC_CDS_PARTS["v0.2"]) == 3
    for value in bsb.BASIC_CDS_PARTS["v0.2"].values():
        assert type(value) == bsb.BasicPart


def test_basic_promoter_parts_v02():
    assert len(bsb.BASIC_PROMOTER_PARTS["v0.2"]) == 60
    for value in bsb.BASIC_PROMOTER_PARTS["v0.2"].values():
        assert type(value) == bsb.BasicPart


def test_basic_promoter_parts_v03():
    assert len(bsb.BASIC_PROMOTER_PARTS["v0.3"]) == 60
    for value in bsb.BASIC_PROMOTER_PARTS["v0.3"].values():
        assert type(value) == bsb.BasicPart


def test_basic_seva_parts_v10():
    assert len(bsb.BASIC_SEVA_PARTS["v1.0"]) == 30
    for value in bsb.BASIC_SEVA_PARTS["v1.0"].values():
        assert type(value) == bsb.BasicPart


def test_linker_384_plate():
    from basicsynbio.cam.default_linker_plate import LINKER_384_PLATE
    from platemap import Plate

    assert type(LINKER_384_PLATE) == Plate
    assert LINKER_384_PLATE.size == 384
    assert LINKER_384_PLATE.deadspace == 20e3
    assert LINKER_384_PLATE.well_volume == 65e3
    for well_id in ("A1", "B1", "A2", "B2"):
        assert "L1-S" in LINKER_384_PLATE[well_id]["composition"].keys()
    for well_id in ("O21", "P21", "O22", "P22"):
        assert "UTR3-S" in LINKER_384_PLATE[well_id]["composition"].keys()
