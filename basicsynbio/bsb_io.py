from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs, IP_SEQREC, IS_SEQREC, seqrec2part
import basicsynbio as bsb
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import requests
import io


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
def import_ice_part(ice_client, ice_token, ice_num, file_type="original", format="genbank"):
    """Returns a BasicPart object using an entry on the JBEI-ICE public registry.

    Args:
        ice_client -- X-ICE-API-Token-Client
        ice_token -- X-ICE-API-Token
        ice_num -- Number of Part ID. For instance, ice_num=17338 for JPUB_017338.
        file_type -- The file type to download e.g. "original", "genbank" or "fasta"."""
    ice_url = f"https://public-registry.jbei.org/rest/file/{ice_num}/sequence/{file_type}"
    try:
        ice_response = requests.get(
            ice_url,
            headers={
                "X-ICE-API-Token-Client": ice_client,
                "X-ICE-API-Token": ice_token,
                "Cache-Control": "no-cache"
            },
            timeout=30
        )
    except TimeoutError:
        print("No response from public-registry.jbei.org after 30 seconds.")
    memory_file = io.StringIO(ice_response.text)
    return seqrec2part(SeqIO.read(memory_file, format))