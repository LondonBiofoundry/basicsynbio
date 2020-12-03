"""Module contains objects for importing and exporting parts and sequences."""

from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs, IP_SEQREC, IS_SEQREC, seqrec2part
import basicsynbio as bsb
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import io
import json
import os
import tempfile
from sbol2 import Document


@add2docs(CommonArgDocs.HANDLE, CommonArgDocs.FORMAT, CommonArgDocs.ADD_I_SEQS)
def import_part(handle, format, add_i_seqs=False):
    """:return: Part object using Bio.SeqIO.read().

    Refer to Biopython documentation for further information on Bio.SeqIO.read().

    :rtype: :py:class:`BasicPart`

    """
    seqrec = SeqIO.read(handle, format)
    return seqrec2part(seqrec, add_i_seqs)


@add2docs(CommonArgDocs.ADD_I_SEQS)
def import_sbol_part(path, add_i_seqs=False):
    """:return: Part object using sbol2.Document.exportToFormat.

    Refer to pysbol2 documentation for further information.

    :rtype: :py:class:`BasicPart` object.
    :param string path: Path to sbol file.

    """
    doc = Document(path)
    fp = tempfile.NamedTemporaryFile(delete=False)
    doc.exportToFormat("GenBank", fp.name)
    seqrec = SeqIO.read(fp.name, "genbank")
    fp.close()
    os.unlink(fp.name)
    return seqrec2part(seqrec, add_i_seqs)


@add2docs(CommonArgDocs.HANDLE, CommonArgDocs.FORMAT, CommonArgDocs.ADD_I_SEQS)
def import_parts(handle, format, add_i_seqs=False):
    """:return: Generator of :py:class:`BasicPart` instances using Bio.SeqIO.parse().

    Refer to Biopython documentation for further information on Bio.SeqIO.parse().

    """
    seqrecs = SeqIO.parse(handle, format)
    yield from (seqrec2part(seqrec, add_i_seqs) for seqrec in seqrecs)


@add2docs(
    CommonArgDocs.HANDLE,
    CommonArgDocs.FORMAT,
)
def export_sequences_to_file(sequences, handle, format="genbank", molecule_type="DNA"):
    """Exports sequences using Bio.SeqIO.write().

    Refer to Biopython documentation for further information on Bio.SeqIO.write().

    :param sequences: objects to export to file handle.
    :type sequences: A single object or iterable of type/s BasicPart, BasicAssembly or Bio.SeqRecord.SeqRecord.
    """
    if type(sequences) in [bsb.BasicPart, bsb.BasicAssembly, SeqRecord]:
        SeqIO.write(_process_basic_object(sequences, molecule_type), handle, format)
    elif hasattr(sequences, "__iter__"):
        sequences = (
            _process_basic_object(basic_object, molecule_type)
            for basic_object in sequences
        )
        SeqIO.write(sequences, handle, format)
    else:
        raise ValueError(
            "sequences was not iterable or of type BasicPart, BasicAssembly or SeqRecord"
        )


def _process_basic_object(basic_object, molecule_type):
    """Converts basic_object into an object that can be processed by Bio.SeqIO."""
    try:
        basic_object = basic_object.return_seqrec()
    except AttributeError:
        pass
    basic_object.annotations["molecule_type"] = molecule_type
    return basic_object
