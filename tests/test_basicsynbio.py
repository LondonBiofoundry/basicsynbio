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
    gfp_orf_feature = feature_from_qualifier(gfp_seqrec, "gene", ["sfGFP"])
    return gfp_orf_feature.extract(gfp_seqrec.seq)

@pytest.fixture
def cmr_p15a_basicpart():
    return basicsynbio.import_basic_part(
        "BASIC_SEVA_37_CmR-p15A.1.gb", "genbank"
    )


@pytest.fixture
def cmr_p15a_backbone():
    from Bio import SeqIO
    cmr_p15a_backbone = SeqIO.read("BASIC_SEVA_37_CmR-p15A.1.gb", "genbank")
    prefix = feature_from_qualifier(cmr_p15a_backbone, "label", ["Prefix"])
    suffix = feature_from_qualifier(cmr_p15a_backbone, "label", ["Suffix"])
    return cmr_p15a_backbone[int(prefix.location.end):] \
        + cmr_p15a_backbone[:int(suffix.location.start)] 


def feature_from_qualifier(seqrec, qualifier_key, qualifier_value):
    """
    Extract the feature from the seqrec that contains the corresponding
    qualifier key/qualifier value pair. Features that

    """
    for feature in seqrec.features:
        if qualifier_key in feature.qualifiers:
            if feature.qualifiers[qualifier_key] == qualifier_value:
                return feature
    raise ValueError(f"{seqrec.id} lacks a feature containing a {qualifier_key}/{qualifier_value} pair in qualifiers")


def test_basic_part(gfp_basicpart, gfp_seqrec):
    """
    Cannot compare features as each is an object
    """
    for key, value in gfp_seqrec.__dict__.items():
        if key != "features":
            assert value == getattr(gfp_basicpart, key)


def test_basic_slice_ip(gfp_basicpart, gfp_orf_seq):
    sliced_part = gfp_basicpart.basic_slice()
    assert sliced_part.seq == gfp_orf_seq


def test_basic_slice_is(cmr_p15a_basicpart, cmr_p15a_backbone):
    sliced_cmr_p15a = cmr_p15a_basicpart.basic_slice()
    assert sliced_cmr_p15a.seq == cmr_p15a_backbone.seq


def test_basic_part_exception(gfp_orf_seq):
    with pytest.raises(basicsynbio_exceptions.PartException):
        basicsynbio.BasicPart(gfp_orf_seq, "sfGFP")

