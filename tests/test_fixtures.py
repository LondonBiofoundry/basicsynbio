import basicsynbio as bsb
import pytest


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


@pytest.fixture
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


@pytest.fixture
def promoter_assemblies_build_more_than_384():
    utr_linkers = [
        "UTR1-RBS1",
        "UTR1-RBS2",
        "UTR1-RBS3",
        "UTR2-RBS1",
        "UTR2-RBS2",
        "UTR2-RBS3",
        "UTR3-RBS1",
        "UTR3-RBS2",
        "UTR3-RBS3",
    ]
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


@pytest.fixture
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
