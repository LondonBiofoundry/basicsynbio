# Uses pytest
from basicsynbio import basicsynbio
import pytest
from Bio import SeqRecord
from Bio import SeqIO


@pytest.fixture
def generic_basic_part():
    return basicsynbio.import_basic_part("generic_basic_part.gb", "genbank")


@pytest.fixture
def generic_seqrec():
    return SeqIO.read("generic_basic_part.gb", "genbank")


def test_basic_part(generic_basic_part, generic_seqrec):
    assert generic_basic_part.seq == generic_seqrec.seq
    assert generic_basic_part.id == generic_seqrec.id


def test_basic_slice_ip_first(generic_basic_part):
    sliced_part = generic_basic_part.basic_slice()
    presliced_part = SeqIO.read("sliced_basic_part_ip.gb", "genbank")
    assert sliced_part.seq == presliced_part.seq

def test_basic_slice_is_first():
    part = basicsynbio.import_basic_part("is_first_part.gb", "genbank")
    sliced_part = part.basic_slice()
    presliced_part = SeqIO.read("sliced_basic_part_is.gb", "genbank")
    assert sliced_part.seq == presliced_part.seq

