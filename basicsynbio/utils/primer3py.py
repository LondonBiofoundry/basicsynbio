"""This module contains objects which make it easier to work with and enhance primer3-py."""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from datetime import datetime

DATE = datetime.now()
OLIGO_ANNOTATIONS = {
    "date": DATE.strftime("%d-") + DATE.strftime("%b").upper() + DATE.strftime("-%Y"),
    "sequence_version": 1,
    "topology": "linear",
    "molecule_type": "DNA",
}


def p3_seqrec(
    p3py_out: dict,
    primer: str,
    seqrec_attrs: dict = None,
    primer_pair: int = 0,
) -> SeqRecord:
    """Return a SeqRecord of a primer3-py primer.

    Args:
        p3py_out: Refer to __init__().
        primer: Either "LEFT" or "RIGHT".
        primer_pair: Refer to __init__().
        seqrec_attrs: Assign SeqRecord attributes. If None, assigns default attributes.

    """
    primer_str = f"PRIMER_{primer}_{str(primer_pair)}"
    seqrec = SeqRecord(
        Seq(p3py_out[primer_str + "_SEQUENCE"]),
        primer_str,
        annotations=OLIGO_ANNOTATIONS,
    )
    if seqrec_attrs:
        for key, value in seqrec_attrs.items():
            setattr(seqrec, key, value)
    return seqrec
