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


def test_return_seqrec(five_part_assembly):
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
