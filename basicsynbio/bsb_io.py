"""Module contains objects for importing and exporting parts and sequences."""

from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs, IP_SEQREC, IS_SEQREC, seqrec2part
import basicsynbio as bsb
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import io
import icebreaker
import json
import os
import tempfile
from sbol2 import Document


@add2docs(
    CommonArgDocs.HANDLE,
    CommonArgDocs.FORMAT,
    CommonArgDocs.ADD_I_SEQS
)
def import_part(handle, format, add_i_seqs=False):
    """:return: Part object using Bio.SeqIO.read().

    Refer to BioPython documentation for further information.

    :rtype: :py:class:`BasicPart`

    """
    seqrec = SeqIO.read(handle, format)
    return seqrec2part(seqrec, add_i_seqs)


@add2docs(
    CommonArgDocs.ADD_I_SEQS
)
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


@add2docs(
    CommonArgDocs.HANDLE,
    CommonArgDocs.FORMAT,
    CommonArgDocs.ADD_I_SEQS
)
def import_parts(handle, format, add_i_seqs=False):
    """:return: Generator of :py:class:`BasicPart` instances using Bio.SeqIO.parse().

    Refer to BioPython documentation for further information.

    """
    seqrecs = SeqIO.parse(handle, format)
    yield from (seqrec2part(seqrec, add_i_seqs) for seqrec in seqrecs)


@add2docs(
    CommonArgDocs.HANDLE,
    CommonArgDocs.FORMAT,
)
def export_sequences_to_file(sequences, handle, format="genbank"):
    """Exports sequences using Bio.SeqIO.write().

    Refer to Bio.SeqIO.write documentation for further information.

    :param sequences: objects to export to file handle. 
    :type sequences: A single object or iterable of type/s :py:class:`BasicPart`, :py:class:`BasicAssembly` or Bio.SeqRecord.SeqRecord.

    """
    if type(sequences) in [bsb.BasicPart, bsb.BasicAssembly, SeqRecord]:
        basic_object = sequences
        try:
            basic_object = basic_object.return_seqrec()
            SeqIO.write(basic_object, handle, format)
        except AttributeError:
            SeqIO.write(basic_object, handle, format)
    else:
        sequences = ((basic_object if not hasattr(basic_object, "return_seqrec")
                      else basic_object.return_seqrec()) for basic_object in sequences)
        SeqIO.write(sequences, handle, format)


def import_ice_parts(
        ice_user_config: dict,
        *ice_nums,
        ice_root="https://public-registry.jbei.org/",
        file_type="original",
        format="genbank"):
    """Returns a BasicPart object using an entry on a JBEI-ICE instance. Uses icebreaker under the hood.

    Refer to icebreaker documentation for further information.

    :param dict ice_user_config: Either {email: password:} or {token: client:}. Note, compared to icebreaker, 'Root' is a separate default argument (ice_root).  
    :param \*ice_nums: Part ID numbers as strings e.g. '17338'. Note, the distinction between Part ID and Part ID number. For instance, number is 17338 for ID=JPUB_017338. 
    :param string ice_root: Root of ice registry.
    :param string file_type: The file type to download e.g. "original", "genbank" or "fasta".
    :param string format: Format of the downloaded file. If :param file_type: == "fasta", :param format: must be correspondingly set to "fasta".
        
    """
    ice_config = ice_user_config
    ice_config["root"] = ice_root
    ice = icebreaker.IceClient(ice_config)
    for ice_num in ice_nums:
        bytes_file = ice.request(
            method="GET",
            endpoint=f"file/{ice_num}/sequence/{file_type}",
            response_type="file"
        )
        memory_file = io.StringIO(bytes_file.decode("utf-8"))
        yield import_part(memory_file, "genbank")