# Uses pytest
from basicsynbio import basicsynbio, basicsynbio_exceptions
import pytest


@pytest.fixture
def gfp_basicpart():
    return basicsynbio.import_basic_part("BASIC_sfGFP_ORF.1.gb", "genbank")


@pytest.fixture
def gfp_seqrec():
    from Bio import SeqIO
    return SeqIO.read("BASIC_sfGFP_ORF.1.gb", "genbank")

@pytest.fixture
def gfp_orf_seq(gfp_seqrec):
    for feature in gfp_seqrec.features:
        if feature.type == "CDS":
            gfp_orf_feature = feature
    return gfp_orf_feature.extract(gfp_seqrec.seq)


def test_basic_part(gfp_basicpart, gfp_seqrec):
    """
    Cannot compare features as each is an object
    """
    for key, value in gfp_seqrec.__dict__.items():
        if key not in ["features"]:
            assert value == getattr(gfp_basicpart, key)


def test_basic_slice_ip_first(gfp_basicpart, gfp_orf_seq):
    sliced_part = gfp_basicpart.basic_slice()
    assert sliced_part.seq == gfp_orf_seq

def test_basic_slice_is_first():
    part = basicsynbio.import_basic_part("is_first_part.gb", "genbank")
    sliced_part = part.basic_slice()

    assert sliced_part.seq == presliced_part.seq


def test_basic_part_exception():
    with pytest.raises(basicsynbio_exceptions.PartException):
        basicsynbio.import_basic_part("gfp_orf", "genbank")

