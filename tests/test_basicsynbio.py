import basicsynbio as bsb
import pytest


class ComparisonException(Exception):
    pass


@pytest.fixture
def gfp_basicpart():
    return bsb.import_part(
        "sequences/genbank_files/misc_BASIC/BASIC_sfGFP_ORF.1.gb", "genbank"
    )


@pytest.fixture
def gfp_seqrec():
    from Bio import SeqIO

    return SeqIO.read(
        "sequences/genbank_files/misc_BASIC/BASIC_sfGFP_ORF.1.gb", "genbank"
    )


@pytest.fixture
def gfp_orf_seq(gfp_seqrec):
    from basicsynbio.utils import feature_from_qualifier

    gfp_orf_feature = feature_from_qualifier(gfp_seqrec, "gene", ["sfGFP"])
    return gfp_orf_feature.extract(gfp_seqrec.seq)


@pytest.fixture
def cmr_p15a_basicpart():
    return bsb.import_part(
        "sequences/genbank_files/previous_versions/BASIC_SEVA_37_CmR-p15A.1.gb",
        "genbank",
    )


@pytest.fixture
def cmr_p15a_backbone():
    from basicsynbio.utils import feature_from_qualifier
    from Bio import SeqIO

    cmr_p15a_backbone = SeqIO.read(
        "sequences/genbank_files/previous_versions/BASIC_SEVA_37_CmR-p15A.1.gb",
        "genbank",
    )
    prefix = feature_from_qualifier(cmr_p15a_backbone, "label", ["Prefix"])
    suffix = feature_from_qualifier(cmr_p15a_backbone, "label", ["Suffix"])
    return (
        cmr_p15a_backbone[int(prefix.location.end) :]
        + cmr_p15a_backbone[: int(suffix.location.start)]
    )


@pytest.fixture
def five_part_assembly_parts(cmr_p15a_basicpart, gfp_basicpart):
    promoter = bsb.import_part(
        "sequences/genbank_files/misc_BASIC/BASIC_L3S2P21_J23105_RiboJ.1.gb", "genbank"
    )
    bfp_basicpart = bsb.import_part(
        "sequences/genbank_files/misc_BASIC/BASIC_mTagBFP2_ORF.1.gb", "genbank"
    )
    rfp_basicpart = bsb.import_part(
        "sequences/genbank_files/misc_BASIC/BASIC_mCherry_ORF.1.gb", "genbank"
    )
    return [cmr_p15a_basicpart, promoter, gfp_basicpart, bfp_basicpart, rfp_basicpart]


@pytest.fixture
def five_part_assembly_linkers():
    return [
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR1-RBS2"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR2-RBS1"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR3-RBS1"],
    ]


@pytest.fixture
def five_part_assembly(five_part_assembly_parts, five_part_assembly_linkers):
    zipped_parts_linkers = zip(five_part_assembly_linkers, five_part_assembly_parts)
    parts_linkers = []
    for part_linker in zipped_parts_linkers:
        parts_linkers += list(part_linker)
    print([part_linker.id for part_linker in parts_linkers])
    print([type(part_linker) for part_linker in parts_linkers])
    return bsb.BasicAssembly("5_part", *parts_linkers)


@pytest.fixture
def gfp_orf_seqrec(gfp_orf_seq):
    from basicsynbio.utils import _easy_seqrec

    return _easy_seqrec(
        str(gfp_orf_seq),
        "sfGFP",
        annotation_type="CDS",
        note=["fluorescent reporter protein"],
        gene=["sfGFP"],
    )


@pytest.fixture
def gfp_orf_basicpart(gfp_orf_seqrec):
    return bsb.seqrec2part(gfp_orf_seqrec, add_i_seqs=True)


@pytest.fixture
def bseva_68_seqrec():
    return bsb.BASIC_SEVA_PARTS["v0.1"]["68"]


@pytest.fixture
def bsai_part_seqrec(gfp_orf_seq):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    bsai_site = Seq("GGTCTC")
    insertion_ind = len(gfp_orf_seq) // 2
    return SeqRecord(
        gfp_orf_seq[:insertion_ind] + bsai_site + gfp_orf_seq[insertion_ind:],
        id="bsai_part",
    )


@pytest.yield_fixture
def promoter_assemblies_build():
    utr_linkers = ["UTR1-RBS1", "UTR1-RBS2", "UTR1-RBS3"]
    promoter_assemblies = []
    for utr_linker in utr_linkers:
        promoter_assemblies += [
            bsb.BasicAssembly(
                f"promoter_construct_{ind}_{utr_linker}",
                bsb.BASIC_SEVA_PARTS["v0.1"]["26"],
                bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
                promoter,
                bsb.BASIC_BIOLEGIO_LINKERS["v0.1"][utr_linker],
                bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
                bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
            )
            for ind, promoter in enumerate(bsb.BASIC_PROMOTER_PARTS["v0.1"].values())
        ]
    return bsb.BasicBuild(*promoter_assemblies)


@pytest.yield_fixture
def all_promoter_assemblies_build():
    promoter_assemblies = []
    promoter_assemblies += [
        bsb.BasicAssembly(
            f"promoter_construct_{ind}",
            bsb.BASIC_SEVA_PARTS["v0.1"]["26"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
            promoter,
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR1-RBS1"],
            bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
        )
        for ind, promoter in enumerate(bsb.BASIC_PROMOTER_PARTS["v0.1"].values())
    ]
    return bsb.BasicBuild(*promoter_assemblies)


@pytest.fixture
def promoter_assemblies_json(promoter_assemblies_build):
    import json

    return json.dumps(promoter_assemblies_build, cls=bsb.BuildEncoder, indent=4)


@pytest.fixture
def gfp_part_final_conc(gfp_basicpart):
    from Bio.SeqUtils import molecular_weight

    return (
        2.5
        * molecular_weight(gfp_basicpart.seq, double_stranded=True, circular=True)
        / 1e6
    )


@pytest.fixture
def small_build_example():
    return bsb.BasicBuild(
        bsb.BasicAssembly(
            "First_Assembly_With_18",
            bsb.BASIC_SEVA_PARTS["v0.1"]["18"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
            bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
        ),
        bsb.BasicAssembly(
            "Second_Assembly_With_26",
            bsb.BASIC_SEVA_PARTS["v0.1"]["26"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
            bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
        ),
    )


def compare_seqrec_instances(seqrec1, seqrec2):
    """
    returns true if seqrec1 has equivalent seqrec2 attributes.
    Ignores seqrec.features as contains SeqFeature objects.

    """
    for key, value in seqrec2.__dict__.items():
        if key != "features":
            if value != getattr(seqrec1, key):
                print(
                    f"seqrec2's '{key}' attribute does not match that obtained from seqrec1."
                )
                return False
    return True


def json_round_trip(class_instance, encoder, decoder):
    """Encodes object as serialised json and returns decoded object.

    Args:
        object -- class instance to be serialised.
        encoder -- json encoder Class.
        decoder -- json decoder function.
    """
    import json

    serialised_object = json.dumps(class_instance, cls=encoder)
    return json.loads(serialised_object, cls=decoder)


def test_basic_part(gfp_basicpart, gfp_seqrec):
    assert compare_seqrec_instances(gfp_basicpart, gfp_seqrec) == True


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


def test_basic_part_final_concentration(gfp_basicpart, gfp_part_final_conc):
    assert gfp_basicpart.concentration(stock=False) == round(gfp_part_final_conc)


def test_basic_part_stock_concentration(gfp_basicpart, gfp_part_final_conc):
    clip_vol = 30
    assert gfp_basicpart.concentration(stock=True) == round(
        gfp_part_final_conc * clip_vol
    )


def test_assembly_error(gfp_basicpart, cmr_p15a_basicpart):
    import basicsynbio.main as bsb_main

    with pytest.raises(bsb_main.AssemblyException):
        bsb.BasicAssembly("test", gfp_basicpart, cmr_p15a_basicpart)


def testreturn_seqrec(five_part_assembly):
    from Bio import SeqIO

    example_assembly = SeqIO.read(
        "sequences/genbank_files/misc_BASIC/five_part_assembly.gb", "genbank"
    )
    assert five_part_assembly.return_seqrec().seq == example_assembly.seq


@pytest.mark.skip(
    reason="Added bsb_io module and removed BasicAssembly.return_file() method"
)
def test_assembly_return_file(five_part_assembly):
    """The BASIC assembly return_file() method is required given all BASIC assemblies might not be BASIC parts."""
    import os

    assembly_seqrec = five_part_assembly.return_seqrec()
    five_part_assembly.return_file("test_five_part_assembly.gb")
    os.remove("test_five_part_assembly.gb")


def test_basic_parts_in_file():
    parts = bsb.import_parts(
        "sequences/genbank_files/misc_BASIC/dnabot_constructs.gb", "genbank"
    )
    print(list(parts)[:5])


def test_add_i_seqs(gfp_orf_basicpart, gfp_orf_seq):
    import basicsynbio.main as bsb_main

    print("length of gfp_basicpart: ", len(gfp_orf_basicpart))
    print(
        "length of correct sequence: ",
        len(bsb_main.IP_STR) + len(gfp_orf_basicpart) + len(bsb_main.IS_STR),
    )
    assert (
        str(gfp_orf_basicpart.seq)
        == bsb_main.IP_STR + str(gfp_orf_seq) + bsb_main.IS_STR
    )
    assert len(gfp_orf_basicpart.features) == 3


def test_return_part(five_part_assembly):
    imported_part = bsb.import_part(
        "sequences/genbank_files/misc_BASIC/five_part_assembly.gb", "genbank"
    )
    api_part = five_part_assembly.return_part()
    assert api_part.seq == imported_part.seq
    assert dir(api_part) == dir(imported_part)


def test_export_to_file(gfp_basicpart, five_part_assembly, gfp_seqrec):
    import os

    try:
        bsb.export_sequences_to_file(gfp_basicpart, "test_export.gb", "genbank")
        print("finished exporting BasicPart")
        bsb.export_sequences_to_file(five_part_assembly, "test_export.gb", "genbank")
        print("finished exporting BasicAssembly")
        bsb.export_sequences_to_file(gfp_seqrec, "test_export.gb", "genbank")
        print("finished exporting SeqRecord")
        bsb.export_sequences_to_file(
            [gfp_basicpart, five_part_assembly, gfp_seqrec], "test_export.gb", "genbank"
        )
        print("finished exporting iterable")
    finally:
        os.remove("test_export.gb")


def test_export_new_part(gfp_orf_seq):
    import os
    from Bio.SeqRecord import SeqRecord

    try:
        seqrec = SeqRecord(gfp_orf_seq, "gfp_orf")
        template = bsb.seqrec2part(seqrec, add_i_seqs=True)
        bsb.export_sequences_to_file(template, "test_export.gb")
        print("finished exporting GenBank")
    finally:
        os.remove("test_export.gb")


def test_export_csv(promoter_assemblies_build):
    import os

    try:
        promoter_assemblies_build.export_csvs("test_build.zip")
        print("finished exporting Assemblies.csv")
        print("finished exporting Clips.csv")
    finally:
        os.remove("test_build.zip")


def test_addargs2docs_decorator():
    from basicsynbio.main import CommonArgDocs
    from basicsynbio.decorators import addargs2docs

    @addargs2docs(CommonArgDocs.ADD_I_SEQS)
    def dummy_func():
        """add_i_seqs:"""
        pass

    print(bsb.seqrec2part.__doc__)
    print(dummy_func.__doc__)
    assert (
        dummy_func.__doc__
        == """add_i_seqs: if True adds flanking BASIC iP and iS sequences. Note, letter_annotations attribute is lost."""
    )


def test_new_part_resuspension(gfp_orf_basicpart):
    from Bio.SeqUtils import molecular_weight

    print(f"length of basicpart: {len(gfp_orf_basicpart.seq)}")
    print(f"estimated MW: {len(gfp_orf_basicpart.seq*660)}")
    print(
        f"biopython MW: {molecular_weight(gfp_orf_basicpart.seq, double_stranded=True)}"
    )
    mass = 750
    vol = bsb.new_part_resuspension(part=gfp_orf_basicpart, mass=mass)
    print(f"Calculated volume of resuspension buffer: {vol}")
    mw = molecular_weight(gfp_orf_basicpart.seq, double_stranded=True)
    print(f"estimated concentration: {mass*1e-9/(vol*1e-6*mw)*1e9}")
    assert 75 == round(mass * 1e-9 / (vol * 1e-6 * mw) * 1e9)


def test_bseva_dict(bseva_68_seqrec):
    print(bsb.BASIC_SEVA_PARTS["v0.1"].keys())
    bseva_68_part = bsb.BASIC_SEVA_PARTS["v0.1"]["68"]
    assert compare_seqrec_instances(bseva_68_part, bseva_68_seqrec) == True


def test_bpromoter_dict():
    from Bio import SeqIO

    bpromoters_handle = "./basicsynbio/parts_linkers/BASIC_promoter_collection.gb"
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


def test_all_feature_values(gfp_orf_basicpart):
    from basicsynbio.utils import all_feature_values

    print(all_feature_values(gfp_orf_basicpart))
    assert all_feature_values(gfp_orf_basicpart) == [
        "BASIC integrated prefix",
        "fluorescent reporter protein",
        "sfGFP",
        "BASIC integrated suffix",
    ]


def test_multiple_integrated_sequences(gfp_orf_seqrec):
    from basicsynbio.main import IP_SEQREC, IS_SEQREC, PartException, seqrec2part

    with pytest.raises(PartException):
        seqrec2part(IP_SEQREC + IP_SEQREC + gfp_orf_seqrec + IS_SEQREC)


def test_BasicAssembly_clips(
    five_part_assembly, five_part_assembly_parts, five_part_assembly_linkers
):
    from basicsynbio.main import ClipReaction

    clips = []
    for ind, part in enumerate(five_part_assembly_parts):
        if ind == len(five_part_assembly_parts) - 1:
            clips.append(
                ClipReaction(
                    prefix=five_part_assembly_linkers[ind],
                    part=part,
                    suffix=five_part_assembly_linkers[0],
                )
            )
        else:
            clips.append(
                ClipReaction(
                    prefix=five_part_assembly_linkers[ind],
                    part=part,
                    suffix=five_part_assembly_linkers[ind + 1],
                )
            )
    for clip in clips:
        assert clip in five_part_assembly.clip_reactions


@pytest.mark.skip(reason="people should realise this is a bad idea!")
def test_assembly_monkey_clips(five_part_assembly):
    five_part_assembly.parts_linkers = (
        bsb.BASIC_SEVA_PARTS["v0.1"]["18"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
        bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
    )
    assert len(five_part_assembly.clip_reactions) == len(
        [
            part
            for part in five_part_assembly.parts_linkers
            if isinstance(part, bsb.BasicPart)
        ]
    )


def test_assembly_exception_same_utr_linker(cmr_p15a_basicpart, gfp_basicpart):
    from basicsynbio.main import AssemblyException

    with pytest.raises(
        AssemblyException, match="BasicAssembly initiated with UTR1-S used 2 times."
    ):
        bsb.BasicAssembly(
            "test",
            cmr_p15a_basicpart,
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR1-RBS1"],
            gfp_basicpart,
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR1-RBS2"],
        )


def test_bsai_site_in_part(bsai_part_seqrec):
    with pytest.raises(
        bsb.main.PartException,
        match=f"{bsai_part_seqrec.id} contains more than two BsaI sites.",
    ):
        bsb.seqrec2part(bsai_part_seqrec, add_i_seqs=True)


def test_build_parts(promoter_assemblies_build):
    parts = [
        promoter_part for promoter_part in bsb.BASIC_PROMOTER_PARTS["v0.1"].values()
    ]
    parts += [bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"], bsb.BASIC_SEVA_PARTS["v0.1"]["26"]]
    print(parts)
    part_ids = [part.id for part in parts]
    for element in promoter_assemblies_build.unique_parts_data.values():
        assert element["part"].id in part_ids


def test_build_linkers(promoter_assemblies_build):
    linkers = ("LMP", "UTR1-RBS1", "UTR1-RBS2", "UTR1-RBS3", "LMS")
    linker_ids = [bsb.BASIC_BIOLEGIO_LINKERS["v0.1"][linker].id for linker in linkers]
    for element in promoter_assemblies_build.unique_linkers_data.values():
        assert element["linker"].id in linker_ids


def test_build_clips_data(promoter_assemblies_build):
    from basicsynbio.main import ClipReaction

    clip_reactions = [
        ClipReaction(
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
            bsb.BASIC_SEVA_PARTS["v0.1"]["26"],
            bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
        )
    ]
    for utr_linker in ("UTR1-RBS1", "UTR1-RBS2", "UTR1-RBS3"):
        clip_reactions.append(
            ClipReaction(
                bsb.BASIC_BIOLEGIO_LINKERS["v0.1"][utr_linker],
                bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
                bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
            )
        )
        clip_reactions += [
            ClipReaction(
                bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
                promoter,
                bsb.BASIC_BIOLEGIO_LINKERS["v0.1"][utr_linker],
            )
            for promoter in bsb.BASIC_PROMOTER_PARTS["v0.1"].values()
        ]
    for clip_reaction in promoter_assemblies_build.clips_data.keys():
        assert clip_reaction in clip_reactions
    assert len(promoter_assemblies_build.clips_data) == len(clip_reactions)


def test_basic_build_indetical_ids(five_part_assembly):
    from basicsynbio.cam import BuildException

    with pytest.raises(
        BuildException,
        match=f"ID '{five_part_assembly.id}' has been assigned to 2 BasicAssembly instance/s. All assemblies of a build should have a unique 'id' attribute.",
    ):
        bsb.BasicBuild(five_part_assembly, five_part_assembly)


def test_unique_parts_in_build_are_unique(promoter_assemblies_build):
    true_unique_parts = []
    for part in promoter_assemblies_build.unique_parts:
        if part not in true_unique_parts:
            true_unique_parts.append(part)
    print(f"true unique part IDs: {[part.id for part in true_unique_parts]}")
    print(
        f"build unique part IDs: {[part.id for part in promoter_assemblies_build.unique_parts]}"
    )
    assert len(true_unique_parts) == len(promoter_assemblies_build.unique_parts)


def test_type_of_unique_parts_is_tuple_of_parts(promoter_assemblies_build):
    assert isinstance(promoter_assemblies_build.unique_parts, tuple)
    assert isinstance(promoter_assemblies_build.unique_parts[0], bsb.BasicPart)


def test_type_of_unique_linkers_is_tuple_of_parts(promoter_assemblies_build):
    assert isinstance(promoter_assemblies_build.unique_linkers, tuple)
    assert isinstance(promoter_assemblies_build.unique_linkers[0], bsb.BasicLinker)


def test_partially_decoded_build(promoter_assemblies_json, promoter_assemblies_build):
    import json

    decoded_build = json.loads(promoter_assemblies_json, cls=bsb.BuildDecoder)
    assert True == isinstance(decoded_build, bsb.BasicBuild)
    assert len(promoter_assemblies_build.basic_assemblies) == len(
        decoded_build.basic_assemblies
    )


def test_decoded_build(promoter_assemblies_build, promoter_assemblies_json):
    import json
    from basicsynbio.cam import seqrecord_hexdigest

    decoded_build = json.loads(promoter_assemblies_json, cls=bsb.BuildDecoder)
    original_parts = (
        part_dict["part"]
        for part_dict in promoter_assemblies_build.unique_parts_data.values()
    )
    decoded_build.update_parts(*original_parts)
    sfgfp_hash = seqrecord_hexdigest(bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"])
    assert (
        compare_seqrec_instances(
            decoded_build.unique_parts_data[sfgfp_hash]["part"],
            bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
        )
        == True
    )


def test_error_raise_basic_slice_less_90():
    from Bio.Seq import Seq

    with pytest.raises(ValueError):
        mypart = bsb.BasicPart(Seq("TCTGGTGGGTCTCTGTCCAAGGCTCGGGAGACCTATCG"), "test")


def test_warning_raise_basic_slice_90_150():
    from Bio.Seq import Seq

    with pytest.warns(UserWarning):
        mypart = bsb.BasicPart(
            Seq(
                "TCTGGTGGGTCTCTGTCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGCTCGGGAGACCTATCG"
            ),
            "test",
        )


def test_echo_instructions_small_build(small_build_example):
    import zipfile
    import os
    import pandas as pd
    import numpy as np
    from pathlib import Path

    echozippath = small_build_example.export_echo_assembly()
    with zipfile.ZipFile(echozippath, "r") as zip_ref:
        try:
            zip_ref.extractall("ECHO_CSVS")
        finally:
            zip_ref.close()
            os.remove(echozippath)
    echo_clips = pd.read_csv(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    echo_water_buffer = pd.read_csv(
        Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv"
    )
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv")
    os.rmdir(Path.cwd() / "ECHO_CSVS")
    expected_clips = [
        ["A1", "A1", 500],
        ["A1", "B1", 500],
        ["B1", "C1", 500],
        ["B1", "B1", 500],
    ]
    expected_water_buffer = [
        ["A1", "A1", 500],
        ["A1", "B1", 3500],
        ["B1", "A1", 500],
        ["B1", "B1", 3500],
    ]
    assert expected_clips == echo_clips.to_numpy().tolist()
    assert expected_water_buffer == echo_water_buffer.to_numpy().tolist()


def test_echo_instructions_small_build_useAllWell_False(small_build_example):
    import zipfile
    import os
    import pandas as pd
    import numpy as np
    from pathlib import Path

    echozippath = small_build_example.export_echo_assembly(useAllWells=False)
    with zipfile.ZipFile(echozippath, "r") as zip_ref:
        try:
            zip_ref.extractall("ECHO_CSVS")
        finally:
            zip_ref.close()
            os.remove(echozippath)
    echo_clips = pd.read_csv(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    echo_water_buffer = pd.read_csv(
        Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv"
    )
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv")
    os.rmdir(Path.cwd() / "ECHO_CSVS")
    expected_clips = [
        ["A1", "A1", 500],
        ["A1", "C1", 500],
        ["B1", "E1", 500],
        ["B1", "C1", 500],
    ]
    expected_water_buffer = [
        ["A1", "A1", 500],
        ["A1", "B1", 3500],
        ["B1", "A1", 500],
        ["B1", "B1", 3500],
    ]
    assert expected_clips == echo_clips.to_numpy().tolist()
    assert expected_water_buffer == echo_water_buffer.to_numpy().tolist()


def test_echo_overflow_wells(all_promoter_assemblies_build):
    import zipfile
    import os
    import pandas as pd
    import numpy as np
    from pathlib import Path

    echozippath = all_promoter_assemblies_build.export_echo_assembly()
    with zipfile.ZipFile(echozippath, "r") as zip_ref:
        try:
            zip_ref.extractall("ECHO_CSVS")
        finally:
            zip_ref.close()
            os.remove(echozippath)
    echo_clips = pd.read_csv(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    echo_water_buffer = pd.read_csv(
        Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv"
    )
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_clips_1.csv")
    os.remove(Path.cwd() / "ECHO_CSVS" / "echo_water_buffer_2.csv")
    os.rmdir(Path.cwd() / "ECHO_CSVS")
    output_as_string = ""
    for item in echo_clips.to_numpy().tolist():
        for subitem in item:
            output_as_string += str(subitem)
    # clips CSV as a string to reduce size cointains well swapping once
    # well has been used 40 times
    expected = "A1A1500A1C1500A1D1500B1A1500B1F1500B1D1500C1A1500C1G1500C1D1500D1A1500D1H1500D1D1500E1A1500E1I1500E1D1500F1A1500F1J1500F1D1500G1A1500G1K1500G1D1500H1A1500H1L1500H1D1500A2A1500A2M1500A2D1500B2A1500B2N1500B2D1500C2A1500C2O1500C2D1500D2A1500D2P1500D2D1500E2A1500E2A2500E2D1500F2A1500F2B2500F2D1500G2A1500G2C2500G2D1500H2A1500H2D2500H2D1500A3A1500A3E2500A3D1500B3A1500B3F2500B3D1500C3A1500C3G2500C3D1500D3A1500D3H2500D3D1500E3A1500E3I2500E3D1500F3A1500F3J2500F3D1500G3A1500G3K2500G3D1500H3A1500H3L2500H3D1500A4A1500A4M2500A4D1500B4A1500B4N2500B4D1500C4A1500C4O2500C4D1500D4A1500D4P2500D4D1500E4A1500E4A3500E4D1500F4A1500F4B3500F4D1500G4A1500G4C3500G4D1500H4A1500H4D3500H4D1500A5A1500A5E3500A5D1500B5A1500B5F3500B5D1500C5A1500C5G3500C5D1500D5A1500D5H3500D5D1500E5A1500E5I3500E5D1500F5A1500F5J3500F5D1500G5A1500G5K3500G5D1500H5A1500H5L3500H5D1500A6B1500A6M3500A6E1500B6B1500B6N3500B6E1500C6B1500C6O3500C6E1500D6B1500D6P3500D6E1500E6B1500E6A4500E6E1500F6B1500F6B4500F6E1500G6B1500G6C4500G6E1500H6B1500H6D4500H6E1500A7B1500A7E4500A7E1500B7B1500B7F4500B7E1500C7B1500C7G4500C7E1500D7B1500D7H4500D7E1500E7B1500E7I4500E7E1500F7B1500F7J4500F7E1500G7B1500G7K4500G7E1500H7B1500H7L4500H7E1500A8B1500A8M4500A8E1500B8B1500B8N4500B8E1500C8B1500C8O4500C8E1500D8B1500D8P4500D8E1500E8B1500E8A5500E8E1500F8B1500F8B5500F8E1500G8B1500G8C5500G8E1500"
    assert expected == output_as_string


def test_echo_instruction_assert_buffer_water_well_errors(small_build_example):
    with pytest.raises(ValueError):
        echozippath = small_build_example.export_echo_assembly(waterWell="D1")
    with pytest.raises(ValueError):
        echozippath = small_build_example.export_echo_assembly(bufferWell="D1")


def test_echo_instruction_assert_buffer_water_well_errors(promoter_assemblies_build):
    # This Build contains ~180 assemblies function should raise errors for builds
    # with more than 96 assemblies
    with pytest.raises(ValueError):
        echozippath = promoter_assemblies_build.export_echo_assembly()


def test_basic_linker_label():
    mylinker = bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"]
    assert "LMP" in mylinker.features[0].qualifiers["label"]


@pytest.mark.slow
def test_import_sbol_part():
    from basicsynbio.cam import seqrecord_hexdigest

    bseva18_from_sbol = next(
        bsb.import_sbol_parts("./sequences/alternative_formats/bseva18.rdf")
    )
    # online converter changes annotations attribute
    bseva18_from_sbol.annotations = bsb.BASIC_SEVA_PARTS["v0.1"]["18"].annotations
    bseva18_from_sbol.id = seqrecord_hexdigest(bseva18_from_sbol)
    assert (
        compare_seqrec_instances(bseva18_from_sbol, bsb.BASIC_SEVA_PARTS["v0.1"]["18"])
        == True
    )
