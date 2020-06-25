from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs, IP_SEQREC, IS_SEQREC, seqrec2part
import basicsynbio as bsb
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


@add2docs(
    8,
    CommonArgDocs.HANDLE,
    CommonArgDocs.FORMAT,
    CommonArgDocs.ADD_I_SEQS
)
def import_part(handle, format, add_i_seqs=False):
    """Return a BASIC part object using Bio.SeqIO.read().
    
    Args:"""
    seqrec = SeqIO.read(handle, format)
    return seqrec2part(seqrec, add_i_seqs)
    

@add2docs(
    8,
    CommonArgDocs.HANDLE,
    CommonArgDocs.FORMAT,
    CommonArgDocs.ADD_I_SEQS
)
def import_parts(handle, format, add_i_seqs=False):
    """Return BASIC parts in a single file using Bio.SeqIO.parse().
    
    Args:"""
    seqrecs = SeqIO.parse(handle, format)
    return [seqrec2part(seqrec, add_i_seqs) for seqrec in seqrecs]


@add2docs(
    8,
    CommonArgDocs.HANDLE,
    CommonArgDocs.FORMAT,
)
def export_to_file(sequences, handle, format="genbank"):
    """Exports sequences using Bio.SeqIO.write().
    
    Args:
        sequences -- objects to export to file handle. BasicPart, BasicAssembly or SeqRecord instances."""
    if type(sequences) in [bsb.BasicPart, bsb.BasicAssembly, SeqRecord]:
        basic_object = sequences
        try:
            basic_object = basic_object.return_seqrec()
            SeqIO.write(basic_object, handle, format)
        except AttributeError:
            SeqIO.write(basic_object, handle, format)
    else:
        sequences = ((basic_object if not hasattr(basic_object, "return_seqrec") else basic_object.return_seqrec()) for basic_object in sequences)
        SeqIO.write(sequences, handle, format)