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
    return [
         cmr_p15a_basicpart,
         promoter,
         gfp_basicpart,
         bfp_basicpart,
         rfp_basicpart
     ]


@pytest.fixture
def five_part_assembly_linkers():
    return [
        bsb.BIOLEGIO_LINKERS["LMS"],
        bsb.BIOLEGIO_LINKERS["LMP"],
        bsb.BIOLEGIO_LINKERS["UTR1-RBS2"], 
        bsb.BIOLEGIO_LINKERS["UTR2-RBS1"],
        bsb.BIOLEGIO_LINKERS["UTR3-RBS1"],
    ]


@pytest.fixture
def five_part_assembly(five_part_assembly_parts, five_part_assembly_linkers):
    zipped_parts_linkers = zip(five_part_assembly_linkers, five_part_assembly_parts)
    parts_linkers = []
    for part_linker in zipped_parts_linkers:
        parts_linkers += list(part_linker)
    print([part_linker.id for part_linker in parts_linkers])
    print([type(part_linker) for part_linker in parts_linkers])
    return bsb.BasicAssembly("five_part_assembly", *parts_linkers)


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
    return bsb.BSEVA_PARTS["68"]


@pytest.fixture
def ice_user_config():
    import os
    ice_client = os.environ.get("JBEI_ICE_CLIENT")    
    ice_token = os.environ.get("JBEI_ICE_TOKEN")
    return {"client": ice_client, "token": ice_token}


@pytest.fixture
def bsai_part_seqrec(gfp_orf_seq):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    bsai_site = Seq("GGTCTC")
    insertion_ind = len(gfp_orf_seq) // 2
    return SeqRecord(
        gfp_orf_seq[:insertion_ind] + bsai_site + gfp_orf_seq[insertion_ind:],
        id="bsai_part"
    )


@pytest.yield_fixture
def promoter_assemblies_build():
    return bsb.BasicBuild(*(bsb.BasicAssembly(
            promoter.id,
            bsb.BSEVA_PARTS["27"],
            bsb.BIOLEGIO_LINKERS["LMP"],
            promoter,
            bsb.BIOLEGIO_LINKERS["UTR1-RBS2"],
            bsb.BCDS_PARTS["sfGFP"],
            bsb.BIOLEGIO_LINKERS["LMS"]
        ) for promoter in bsb.BPROMOTER_PARTS.values()))


@pytest.fixture
def promoter_assemblies_json(promoter_assemblies_build):
    import json
    return json.dumps(promoter_assemblies_build, cls=bsb.BuildEncoder, indent=4)


def compare_seqrec_instances(seqrec1, seqrec2):
    """
    returns true if seqrec1 has equivalent seqrec2 attributes.
    Ignores seqrec.features as contains SeqFeature objects.

    """
    for key, value in seqrec2.__dict__.items():
        if key != "features":
            if value != getattr(seqrec1, key):
                print(f"seqrec2's '{key}' attribute does not match that obtained from seqrec1.")
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


@pytest.mark.skip(reason="Added bsb_io module and removed BasicAssembly.return_file() method")
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
    print(list(parts)[:5])


def test_add_i_seqs(gfp_orf_basicpart, gfp_orf_seq):
    import basicsynbio.main as bsb_main
    print("length of gfp_basicpart: ", len(gfp_orf_basicpart))
    print("length of correct sequence: ", len(bsb_main.IP_STR) + len(gfp_orf_basicpart) + len(bsb_main.IS_STR))
    assert str(gfp_orf_basicpart.seq) == bsb_main.IP_STR + str(gfp_orf_seq) + bsb_main.IS_STR
    assert len(gfp_orf_basicpart.features) == 3


def test_return_part(five_part_assembly):
    imported_part = bsb.import_part("sequences/genbank_files/misc_BASIC/five_part_assembly.gb", "genbank")
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
        bsb.export_sequences_to_file([gfp_basicpart, five_part_assembly, gfp_seqrec], "test_export.gb", "genbank")
        print("finished exporting iterable")
    finally:
        os.remove("test_export.gb")


def test_add2docs_decorator():
    from basicsynbio.main import CommonArgDocs
    core_doc = """Convert SeqRecord to :py:class:`BasicPart`.

    Relevant attributes are maintained.

    :param seqrec: SeqRecord to be converted to :py:class:`BasicPart` subclass.
    :type seqrec: Bio.SeqRecord.SeqRecord
    """
    print(bsb.seqrec2part.__doc__)
    assert bsb.seqrec2part.__doc__ == core_doc + "\n" + " "*4 + CommonArgDocs.ADD_I_SEQS


def test_new_part_resuspension(gfp_orf_basicpart):
    from Bio.SeqUtils import molecular_weight
    print(f"length of basicpart: {len(gfp_orf_basicpart.seq)}")
    print(f"estimated MW: {len(gfp_orf_basicpart.seq*660)}")
    print(f"biopython MW: {molecular_weight(gfp_orf_basicpart.seq, double_stranded=True)}")
    mass = 750
    vol = bsb.new_part_resuspension(part=gfp_orf_basicpart, mass=mass)
    mw = molecular_weight(gfp_orf_basicpart.seq, double_stranded=True)
    print(f"estimated concentration: {mass*1e-9/(vol*1e-6*mw)*1e9}")
    assert 75 == round(mass*1e-9/(vol*1e-6*mw)*1e9)


@pytest.mark.slow
def test_import_ice_parts(bseva_68_seqrec, ice_user_config):
    ice_nums= ["17337"]
    print(f"ice_user_config before import parts: {ice_user_config}")
    ice_parts = bsb.import_ice_parts(ice_user_config, *ice_nums)
    assert compare_seqrec_instances(next(ice_parts), bseva_68_seqrec) == True


@pytest.mark.slow
def test_import_all_ice_parts(ice_user_config):
    ice_nums = (value for value in bsb.BSEVA_ICE_DICT.values())
    ice_parts = bsb.import_ice_parts(ice_user_config, *ice_nums)
    assert len(list(ice_parts)) == len(bsb.BSEVA_ICE_DICT)


def test_bseva_dict(bseva_68_seqrec):
    print(bsb.BSEVA_PARTS.keys())
    bseva_68_part = bsb.BSEVA_PARTS["68"]
    assert compare_seqrec_instances(bseva_68_part, bseva_68_seqrec) == True


def test_bpromoter_dict():
    from Bio import SeqIO
    bpromoters_handle = "./basicsynbio/parts_linkers/BASIC_promoter_collection.gb"
    bpromoter_seqrecs = SeqIO.parse(bpromoters_handle, "genbank")
    for seqrec in bpromoter_seqrecs:
        assert compare_seqrec_instances(bsb.BPROMOTER_PARTS[seqrec.id], seqrec) == True


def test_bcds_dict():
    from Bio import SeqIO
    bcds_handle = "./basicsynbio/parts_linkers/BASIC_CDS_collection.gb"
    bcds_seqrecs = SeqIO.parse(bcds_handle, "genbank")
    for seqrec in bcds_seqrecs:
        assert compare_seqrec_instances(bsb.BCDS_PARTS[seqrec.id], seqrec) == True

    
def test_all_feature_values(gfp_orf_basicpart):
    from basicsynbio.utils import all_feature_values
    print(all_feature_values(gfp_orf_basicpart))
    assert all_feature_values(gfp_orf_basicpart) == ["BASIC integrated prefix", "fluorescent reporter protein", "sfGFP", "BASIC integrated suffix"]


def test_multiple_integrated_sequences(gfp_orf_seqrec):
    from basicsynbio.main import IP_SEQREC, IS_SEQREC, PartException, seqrec2part
    with pytest.raises(PartException):
        seqrec2part(IP_SEQREC + IP_SEQREC + gfp_orf_seqrec + IS_SEQREC)


def test_BasicAssembly_clips(five_part_assembly, five_part_assembly_parts, five_part_assembly_linkers):
    from basicsynbio.main import ClipReaction
    clips = []
    for ind, part in enumerate(five_part_assembly_parts):
        if ind == len(five_part_assembly_parts) - 1:
            clips.append(
                ClipReaction(
                    prefix=five_part_assembly_linkers[ind],
                    part=part,
                    suffix=five_part_assembly_linkers[0]
                )
            )
        else:
            clips.append(
                ClipReaction(
                    prefix=five_part_assembly_linkers[ind],
                    part=part,
                    suffix=five_part_assembly_linkers[ind + 1]
                )
            )
    for clip in clips:
        assert clip in five_part_assembly.clip_reactions
    

@pytest.mark.skip(reason="people should realise this is a bad idea!")
def test_assembly_monkey_clips(five_part_assembly):
    five_part_assembly.parts_linkers = (bsb.BSEVA_PARTS["18"], bsb.BIOLEGIO_LINKERS["LMP"], bsb.BCDS_PARTS["sfGFP"], bsb.BIOLEGIO_LINKERS["LMS"])
    assert len(five_part_assembly.clip_reactions) == len([part for part in five_part_assembly.parts_linkers if isinstance(part, bsb.BasicPart)])


def test_assembly_exception_same_utr_linker(cmr_p15a_basicpart, gfp_basicpart):
    from basicsynbio.main import AssemblyException
    with pytest.raises(AssemblyException, match="BasicAssembly initiated with UTR1-S used 2 times."):
        bsb.BasicAssembly(
            "test",
            cmr_p15a_basicpart,
            bsb.BIOLEGIO_LINKERS["UTR1-RBS1"],
            gfp_basicpart,
            bsb.BIOLEGIO_LINKERS["UTR1-RBS2"]
        )


def test_bsai_site_in_part(bsai_part_seqrec):
    with pytest.raises(bsb.main.PartException, match=f"{bsai_part_seqrec.id} contains more than two BsaI sites."):
        bsb.seqrec2part(bsai_part_seqrec, add_i_seqs=True)


def test_build_parts(promoter_assemblies_build):
    parts = [promoter_part for promoter_part in bsb.BPROMOTER_PARTS.values()]
    parts += [bsb.BCDS_PARTS["sfGFP"], bsb.BSEVA_PARTS["27"]]
    print(parts)
    part_ids = [part.id for part in parts]
    for element in promoter_assemblies_build.unique_parts.values():
        assert element["part"].id in part_ids


def test_build_linkers(promoter_assemblies_build):
    linkers = ("LMP", "UTR1-RBS2", "LMS")
    linker_ids = [bsb.BIOLEGIO_LINKERS[linker].id for linker in linkers]
    for element in promoter_assemblies_build.unique_linkers.values():
        assert element["linker"].id in linker_ids


def test_build_clips_data(promoter_assemblies_build):
    from basicsynbio.main import ClipReaction
    clip_reactions = [
        ClipReaction(bsb.BIOLEGIO_LINKERS["LMS"], bsb.BSEVA_PARTS["27"], bsb.BIOLEGIO_LINKERS["LMP"]),
        ClipReaction(bsb.BIOLEGIO_LINKERS["UTR1-RBS2"], bsb.BCDS_PARTS["sfGFP"], bsb.BIOLEGIO_LINKERS["LMS"]),
    ]
    clip_reactions += [ClipReaction(bsb.BIOLEGIO_LINKERS["LMP"], promoter, bsb.BIOLEGIO_LINKERS["UTR1-RBS2"]) for promoter in bsb.BPROMOTER_PARTS.values()]
    for clip_reaction in promoter_assemblies_build.clips_data.keys():
        assert clip_reaction in clip_reactions
    assert len(promoter_assemblies_build.clips_data) == len(clip_reactions)


def test_basic_build_indetical_ids(five_part_assembly):
    from basicsynbio.cam import BuildException
    with pytest.raises(
        BuildException,
        match=f"ID '{five_part_assembly.id}' has been assigned to 2 BasicAssembly instance/s. All assemblies of a build should have a unique 'id' attribute."
    ):
        bsb.BasicBuild(
            five_part_assembly,
            five_part_assembly
        )


def test_unique_parts_in_build_are_unique(promoter_assemblies_build):
    true_unique_linkers = []
    for linker in promoter_assemblies_build.unique_linkers.values():
        if linker not in true_unique_linkers:
            true_unique_linkers.append(linker)
    assert len(true_unique_linkers) == len(promoter_assemblies_build.unique_linkers)


def test_partially_decoded_build(promoter_assemblies_json, promoter_assemblies_build):
    import json
    decoded_build = json.loads(promoter_assemblies_json, cls=bsb.BuildDecoder)
    assert True == isinstance(decoded_build, bsb.BasicBuild)
    assert len(promoter_assemblies_build.basic_assemblies) == len(decoded_build.basic_assemblies)


def test_decoded_build(promoter_assemblies_build, promoter_assemblies_json):
    import json
    from basicsynbio.cam import _seqrecord_hexdigest
    decoded_build = json.loads(promoter_assemblies_json, cls=bsb.BuildDecoder)
    original_parts = (part_dict["part"] for part_dict in promoter_assemblies_build.unique_parts.values())
    decoded_build.update_parts(*original_parts)
    sfgfp_hash = _seqrecord_hexdigest(bsb.BCDS_PARTS["sfGFP"])
    assert compare_seqrec_instances(decoded_build.unique_parts[sfgfp_hash]["part"], bsb.BCDS_PARTS["sfGFP"]) == True


@pytest.mark.slow
def test_import_sbol_part():
    bseva18_from_sbol = bsb.import_sbol_part("./sequences/alternative_formats/bseva18.rdf")#
    # online converter changes annotations attribute
    bseva18_from_sbol.annotations = bsb.BSEVA_PARTS["18"].annotations
    assert compare_seqrec_instances(bseva18_from_sbol, bsb.BSEVA_PARTS["18"]) == True
