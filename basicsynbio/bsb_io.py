from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs, IP_SEQREC, IS_SEQREC, seqrec2part
import basicsynbio as bsb
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import requests
import io
import icebreaker


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


@add2docs(
    8,
    CommonArgDocs.FORMAT,
)
def import_ice_parts(
    ice_user_config: dict,
    *ice_nums,
    ice_root="https://public-registry.jbei.org/",
    file_type="original",
    format="genbank"):
    """Returns a BasicPart object using an entry on the JBEI-ICE public registry.

    Args:
        ice_user_config -- Either {email: password:} or {token: client:}. Refer to icebreaker documentation. Root is a separate default argument (ice_root).
        *ice_nums -- Part ID numbers e.g. 17338, ... Note, the distinction between Part ID and Part ID number. For instance, number is 17338 for ID=JPUB_017338.
        ice_root -- Root of ice registry. Refer to icebreaker documentation.
        file_type -- The file type to download e.g. "original", "genbank" or "fasta"."""
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
