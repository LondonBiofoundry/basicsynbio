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
