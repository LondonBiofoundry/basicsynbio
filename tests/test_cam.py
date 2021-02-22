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


def test_export_csv(promoter_assemblies_build):
    import os

    try:
        promoter_assemblies_build.export_csvs("test_build.zip")
        print("finished exporting Assemblies.csv")
        print("finished exporting Clips.csv")
    finally:
        os.remove("test_build.zip")


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
