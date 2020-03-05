import basicsynbio as bsb
import basicsynbio.main as bsb_main
import pytest


class ComparisonException(Exception):
    pass


@pytest.fixture
def gfp_basicpart():
    return bsb.import_part("genbank_files/BASIC_sfGFP_ORF.1.gb", "genbank")


@pytest.fixture
def gfp_seqrec():
    from Bio import SeqIO
    return SeqIO.read("genbank_files/BASIC_sfGFP_ORF.1.gb", "genbank")


@pytest.fixture
def gfp_orf_seq(gfp_seqrec):
    gfp_orf_feature = bsb.feature_from_qualifier(gfp_seqrec, "gene", ["sfGFP"])
    return gfp_orf_feature.extract(gfp_seqrec.seq)


@pytest.fixture
def cmr_p15a_basicpart():
    return bsb.import_part(
        "genbank_files/BASIC_SEVA_37_CmR-p15A.1.gb", "genbank"
    )


@pytest.fixture
def cmr_p15a_backbone():
    from Bio import SeqIO
    cmr_p15a_backbone = SeqIO.read("genbank_files/BASIC_SEVA_37_CmR-p15A.1.gb", "genbank")
    prefix = bsb.feature_from_qualifier(cmr_p15a_backbone, "label", ["Prefix"])
    suffix = bsb.feature_from_qualifier(cmr_p15a_backbone, "label", ["Suffix"])
    return cmr_p15a_backbone[int(prefix.location.end):] \
        + cmr_p15a_backbone[:int(suffix.location.start)]


@pytest.fixture
def five_part_assembly(cmr_p15a_basicpart, gfp_basicpart):
     promoter = bsb.import_part(
         "genbank_files/BASIC_L3S2P21_J23105_RiboJ.1.gb", "genbank"
     )
     bfp_basicpart = bsb.import_part(
         "genbank_files/BASIC_mTagBFP2_ORF.1.gb", "genbank"
     )
     rfp_basicpart = bsb.import_part(
         "genbank_files/BASIC_mCherry_ORF.1.gb", "genbank"
     )
     return bsb.BasicAssembly(
         bsb.biolegio_dict["LMS"], cmr_p15a_basicpart, bsb.biolegio_dict["LMP"], \
             promoter, bsb.biolegio_dict["UTR1-RBS2"], gfp_basicpart, \
                 bsb.biolegio_dict["UTR2-RBS1"], bfp_basicpart, \
                     bsb.biolegio_dict["UTR3-RBS1"], rfp_basicpart 
     )


def compare_basicpart_seqrec(basicpart, seqrec):
    """
    returns true if basicpart has equivalent seqrec attributes.
    Ignores seqrec.features as contains objects.

    """
    try:
        for key, value in seqrec.__dict__.items():
            if key != "features":
                if value != getattr(basicpart, key):
                    raise ComparisonException
    except ComparisonException:
        print(f"{basicpart.id} != {seqrec.id}")
        return False
    else:
        return True


def test_basic_part(gfp_basicpart, gfp_seqrec):
    assert compare_basicpart_seqrec(gfp_basicpart, gfp_seqrec) == True


def test_basic_slice_ip(gfp_basicpart, gfp_orf_seq):
    sliced_part = gfp_basicpart.basic_slice()
    assert sliced_part.seq == gfp_orf_seq


def test_basic_slice_is_and_features(cmr_p15a_basicpart, cmr_p15a_backbone):
    sliced_cmr_p15a = cmr_p15a_basicpart.basic_slice()
    print(
        f"sliced backbone features: {sliced_cmr_p15a.features}\nextracted backbone features: {cmr_p15a_backbone.features}"
    )
    assert sliced_cmr_p15a.seq == cmr_p15a_backbone.seq
    assert len(sliced_cmr_p15a.features) == len(cmr_p15a_backbone.features)


def test_basic_part_exception(gfp_orf_seq):
    with pytest.raises(bsb_main.PartException):
        bsb.BasicPart(gfp_orf_seq, "sfGFP")


def test_assembly_error(gfp_basicpart, cmr_p15a_basicpart):
    with pytest.raises(bsb_main.AssemblyException):
        bsb.BasicAssembly(gfp_basicpart, cmr_p15a_basicpart)


def test_return_seqrec(five_part_assembly):
    from Bio import SeqIO
    example_assembly = SeqIO.read(
        "genbank_files/five_part_assembly.gb", "genbank"
    )
    assert five_part_assembly._return_seqrec().seq == example_assembly.seq


def test_assembly_return_file(five_part_assembly):
    import os
    assembly_seqrec = five_part_assembly._return_seqrec()
    print(assembly_seqrec.seq.alphabet)
    five_part_assembly.return_file("test_five_part_assembly.gb")
    os.remove("test_five_part_assembly.gb")


def test_basic_parts_in_file():
    import os
    parts = bsb.import_parts(
        "genbank_files/dnabot_constructs.gb", "genbank"
    )
    print(parts[:5])