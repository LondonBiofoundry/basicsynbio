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

    bcds_handle = "./basicsynbio/parts_linkers/BASIC_CDS_collection.gb"
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


def test_basic_promoter_parts_v02():
    assert len(bsb.BASIC_PROMOTER_PARTS["v0.2"]) == 60
