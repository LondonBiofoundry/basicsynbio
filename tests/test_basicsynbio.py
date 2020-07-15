import basicsynbio as bsb
import pytest


class ComparisonException(Exception):
    pass


@pytest.fixture
def gfp_basicpart():
    return bsb.import_part("sequences/genbank_files/misc_BASIC/BASIC_sfGFP_ORF.1.gb", "genbank")


@pytest.fixture
def gfp_seqrec():
    from Bio import SeqIO
    return SeqIO.read("sequences/genbank_files/misc_BASIC/BASIC_sfGFP_ORF.1.gb", "genbank")


@pytest.fixture
def gfp_orf_seq(gfp_seqrec):
    from basicsynbio.utils import feature_from_qualifier
    gfp_orf_feature = feature_from_qualifier(gfp_seqrec, "gene", ["sfGFP"])
    return gfp_orf_feature.extract(gfp_seqrec.seq)


@pytest.fixture
def cmr_p15a_basicpart():
    return bsb.import_part(
        "sequences/genbank_files/previous_versions/BASIC_SEVA_37_CmR-p15A.1.gb", "genbank"
    )


@pytest.fixture
def cmr_p15a_backbone():
    from basicsynbio.utils import feature_from_qualifier
    from Bio import SeqIO
    cmr_p15a_backbone = SeqIO.read("sequences/genbank_files/previous_versions/BASIC_SEVA_37_CmR-p15A.1.gb", "genbank")
    prefix = feature_from_qualifier(cmr_p15a_backbone, "label", ["Prefix"])
    suffix = feature_from_qualifier(cmr_p15a_backbone, "label", ["Suffix"])
    return cmr_p15a_backbone[int(prefix.location.end):] \
        + cmr_p15a_backbone[:int(suffix.location.start)]


@pytest.fixture
def five_part_assembly(cmr_p15a_basicpart, gfp_basicpart):
     promoter = bsb.import_part(
         "sequences/genbank_files/misc_BASIC/BASIC_L3S2P21_J23105_RiboJ.1.gb", "genbank"
     )
     bfp_basicpart = bsb.import_part(
         "sequences/genbank_files/misc_BASIC/BASIC_mTagBFP2_ORF.1.gb", "genbank"
     )
     rfp_basicpart = bsb.import_part(
         "sequences/genbank_files/misc_BASIC/BASIC_mCherry_ORF.1.gb", "genbank"
     )
     return bsb.BasicAssembly(
         bsb.biolegio_dict["LMS"], cmr_p15a_basicpart, bsb.biolegio_dict["LMP"], \
             promoter, bsb.biolegio_dict["UTR1-RBS2"], gfp_basicpart, \
                 bsb.biolegio_dict["UTR2-RBS1"], bfp_basicpart, \
                     bsb.biolegio_dict["UTR3-RBS1"], rfp_basicpart 
     )


@pytest.fixture
def gfp_orf_basicpart(gfp_orf_seq):
    from basicsynbio.utils import _easy_seqrec
    gfp_orf_seqrec = _easy_seqrec(
            str(gfp_orf_seq),
            "sfGFP",
            annotation_type="CDS",
            note=["fluorescent reporter protein"],
            gene=["sfGFP"],
        )
    return bsb.seqrec2part(gfp_orf_seqrec, add_i_seqs=True)


@pytest.fixture
def bseva_68_seqrec():
    from Bio import SeqIO
    from pathlib import Path
    bseva_dir = Path.cwd() / "sequences" / "genbank_files" / "BASIC_SEVA_collection"
    return SeqIO.read(bseva_dir / "BASIC_SEVA_68.gb", "genbank")


def compare_basicpart_seqrec(basicpart, seqrec):
    """
    returns true if basicpart has equivalent seqrec attributes.
    Ignores seqrec.features as contains SeqFeature objects.

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
    import basicsynbio.main as bsb_main
    with pytest.raises(bsb_main.PartException):
        bsb.BasicPart(gfp_orf_seq, "sfGFP")


def test_assembly_error(gfp_basicpart, cmr_p15a_basicpart):
    import basicsynbio.main as bsb_main
    with pytest.raises(bsb_main.AssemblyException):
        bsb.BasicAssembly(gfp_basicpart, cmr_p15a_basicpart)


def testreturn_seqrec(five_part_assembly):
    from Bio import SeqIO
    example_assembly = SeqIO.read(
        "sequences/genbank_files/misc_BASIC/five_part_assembly.gb", "genbank"
    )
    assert five_part_assembly.return_seqrec().seq == example_assembly.seq


@pytest.mark.skip(reason="Added io module and removed BasicAssembly.return_file() method")
def test_assembly_return_file(five_part_assembly):
    """The BASIC assembly return_file() method is required given all BASIC assemblies might not be BASIC parts."""
    import os
    assembly_seqrec = five_part_assembly.return_seqrec()
    print(assembly_seqrec.seq.alphabet)
    five_part_assembly.return_file("test_five_part_assembly.gb")
    os.remove("test_five_part_assembly.gb")


def test_basic_parts_in_file():
    parts = bsb.import_parts(
        "sequences/genbank_files/misc_BASIC/dnabot_constructs.gb", "genbank"
    )
    print(parts[:5])


def test_add_i_seqs(gfp_orf_basicpart, gfp_orf_seq):
    import basicsynbio.main as bsb_main
    print("length of gfp_basicpart: ", len(gfp_orf_basicpart))
    print("length of correct sequence: ", len(bsb_main.IP_STR) + len(gfp_orf_basicpart) + len(bsb_main.IS_STR))
    assert str(gfp_orf_basicpart.seq) == bsb_main.IP_STR + str(gfp_orf_seq) + bsb_main.IS_STR
    assert len(gfp_orf_basicpart.features) == 3


def test_return_part(five_part_assembly):
    imported_part = bsb.import_part("sequences/genbank_files/misc_BASIC/five_part_assembly.gb", "genbank")
    api_part = five_part_assembly.return_part("five part assembly")
    assert api_part.seq == imported_part.seq
    assert dir(api_part) == dir(imported_part)


def test_export_to_file(gfp_basicpart, five_part_assembly, gfp_seqrec):
    import os
    try:
        bsb.export_to_file(gfp_basicpart, "test_export.gb", "genbank")
        print("finished exporting BasicPart")
        bsb.export_to_file(five_part_assembly, "test_export.gb", "genbank")
        print("finished exporting BasicAssembly")
        bsb.export_to_file(gfp_seqrec, "test_export.gb", "genbank")
        print("finished exporting SeqRecord")
        bsb.export_to_file([gfp_basicpart, five_part_assembly, gfp_seqrec], "test_export.gb", "genbank")
        print("finished exporting iterable")
    finally:
        os.remove("test_export.gb")


def test_add2docs_decorator():
    from basicsynbio.main import CommonArgDocs
    core_doc = "Convert a Bio.SeqRecord to a BasicPart, relevant attributes are maintained.\n\n    Args:\n    "
    print(bsb.seqrec2part.__doc__)
    assert bsb.seqrec2part.__doc__ == core_doc + CommonArgDocs.ADD_I_SEQS


def test_new_part_resuspension(gfp_orf_basicpart):
    from basicsynbio.utils import new_part_resuspension
    from Bio.SeqUtils import molecular_weight
    print(f"length of basicpart: {len(gfp_orf_basicpart.seq)}")
    print(f"estimated MW: {len(gfp_orf_basicpart.seq*660)}")
    print(f"biopython MW: {molecular_weight(gfp_orf_basicpart.seq, double_stranded=True)}")
    mass = 750
    vol = new_part_resuspension(part=gfp_orf_basicpart, mass=mass)
    mw = molecular_weight(gfp_orf_basicpart.seq, double_stranded=True)
    print(f"estimated concentration: {mass*1e-9/(vol*1e-6*mw)*1e9}")
    assert 75 == round(mass*1e-9/(vol*1e-6*mw)*1e9)


def test_import_ice_part(bseva_68_seqrec):
    import os
    bseva_68_ice_number = "17337"
    ice_client = os.environ.get("JBEI_ICE_CLIENT")    
    ice_token = os.environ.get("JBEI_ICE_TOKEN")
    bseva_68_part = bsb.import_ice_part(ice_client, ice_token, bseva_68_ice_number)
    assert compare_basicpart_seqrec(bseva_68_part, bseva_68_seqrec) == True