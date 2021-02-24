import basicsynbio as bsb
import pytest

from .test_fixtures import (
    gfp_basicpart,
    gfp_seqrec,
    gfp_orf_seq,
    cmr_p15a_basicpart,
    cmr_p15a_backbone,
    five_part_assembly_parts,
    five_part_assembly_linkers,
    five_part_assembly,
    gfp_orf_seqrec,
    gfp_orf_basicpart,
    bseva_68_seqrec,
    bsai_part_seqrec,
    promoter_assemblies_build,
    promoter_assemblies_json,
    gfp_part_final_conc,
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
