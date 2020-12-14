"""Module contains objects for importing and exporting parts and sequences."""

from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs, IP_SEQREC, IS_SEQREC, seqrec2part
import basicsynbio as bsb
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import io
import icebreaker
import json
import os
import tempfile
from typing import Union, Iterable, Iterator, Generator
from sbol2 import Document


@add2docs(CommonArgDocs.HANDLE, CommonArgDocs.FORMAT, CommonArgDocs.ADD_I_SEQS)
def import_part(handle: str, format: str, add_i_seqs: bool =False) -> bsb.BasicPart:
    """Imports a Part object using Bio.SeqIO.read().

    Note:
        Refer to Biopython documentation for further information on Bio.SeqIO.read().

    Args:
        handle: handle of file to be parsed
        format: format of handle file could be 'fasta', 'genbank'...
        add_i_seqs (optional): if True adds flanking BASIC iP and iS sequences.
            Note, letter_annotations attribute is lost."

    Returns:
        BasicPart: a Part object using Bio.SeqIO.read()
    """
    seqrec = SeqIO.read(handle, format)
    return seqrec2part(seqrec, add_i_seqs)


@add2docs(CommonArgDocs.ADD_I_SEQS)
def import_sbol_part(path: str, add_i_seqs=False) -> bsb.BasicPart:
    """Imports a BasicPart object using sbol2.Document.exportToFormat.

    Note:
        Refer to Biopython documentation for further information on Bio.SeqIO.read().
        Refer to pysbol2 documentation for further information.

    Args:
        path: path to SBOL file.
        add_i_seqs (optional): if True adds flanking BASIC iP and iS sequences.
            Note, letter_annotations attribute is lost."

    Returns:
        BasicPart: a Part object using Bio.SeqIO.read()
    """
    doc = Document(path)
    fp = tempfile.NamedTemporaryFile(delete=False)
    doc.exportToFormat("GenBank", fp.name)
    seqrec = SeqIO.read(fp.name, "genbank")
    fp.close()
    os.unlink(fp.name)
    return seqrec2part(seqrec, add_i_seqs)


@add2docs(CommonArgDocs.HANDLE, CommonArgDocs.FORMAT, CommonArgDocs.ADD_I_SEQS)
def import_parts(handle: str , format: str, add_i_seqs=False) -> Iterable[bsb.BasicPart]:
    """Imports a Generator of BasicPart objects using Bio.SeqIO.parse().

    Note:
        Refer to Biopython documentation for further information on Bio.SeqIO.parse().

    Args:
        handle: handle of file to be parsed
        format: format of handle file could be 'fasta', 'genbank'...
        add_i_seqs (optional): if True adds flanking BASIC iP and iS sequences.
            Note, letter_annotations attribute is lost."

    Yields:
        BasicPart: all BasicPart objects within the file.
    """
    seqrecs = SeqIO.parse(handle, format)
    yield from (seqrec2part(seqrec, add_i_seqs) for seqrec in seqrecs)


@add2docs(
    CommonArgDocs.HANDLE,
    CommonArgDocs.FORMAT,
)
def export_sequences_to_file(sequences: Iterable[Union[SeqRecord, bsb.BasicPart, bsb.BasicAssembly]], handle: str, format: str ="genbank", molecule_type: str ="DNA") -> None:
    """Exports sequences to file using Bio.SeqIO.write().

    Note:
        Refer to Biopython documentation for further information on Bio.SeqIO.write().

    Args:
        sequences: Sequences to export.
        handle: File name to write to.
        format (optional): Format of handle file could be 'fasta', 'genbank'.
            Defaults to 'genbank'.
        molecule_type (optional): Type of molecule within sequences.
            defaults to 'DNA'.

    Raises:
        ValueError: sequences was not of correct type.

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
        raise TypeError(
            "sequences was not iterable or of type BasicPart, BasicAssembly or SeqRecord"
        )


def _process_basic_object(basic_object, molecule_type):
    """Converts basic_object into an object that can be processed by Bio.SeqIO.

    """
    try:
        basic_object = basic_object.return_seqrec()
    except AttributeError:
        pass
    basic_object.annotations["molecule_type"] = molecule_type
    return basic_object


def import_ice_parts(
    ice_user_config: dict,
    *ice_nums: str,
    ice_root: str ="https://public-registry.jbei.org/",
    file_type: str ="original",
    format: str ="genbank",
) -> Iterator[bsb.BasicPart]:
    """Imports JBEI-ICE instances as BasicPart objects.

    Note:
        Uses icebreaker under the hood. Refer to `icebreaker documentation
            <https://edinburgh-genome-foundry.github.io/icebreaker/>` for 
            further information.
        Note compared to icebreaker, 'Root' is a separate default argument (ice_root).

    Args:
        ice_user_config : Either {email: password:} or {token: client:}
            Note compared to icebreaker, 'Root' is a separate default argument 
            (ice_root).
        *ice_nums: Part ID numbers as strings e.g. "17338". Note, the distinction
            between Part ID and Part ID number. For instance, number is 17338
            for ID=JPUB_017338.
        ice_root (optional): Root of ice registry.
        file_type (optional): The file type to download e.g. "original", "genbank" or "fasta".
        format (optional): Format of the downloaded file. If file_type == "fasta", format must also be "fasta".

    """
    ice_config = ice_user_config
    ice_config["root"] = ice_root
    ice = icebreaker.IceClient(ice_config)
    for ice_num in ice_nums:
        bytes_file = ice.request(
            method="GET",
            endpoint=f"file/{ice_num}/sequence/{file_type}",
            response_type="file",
        )
        memory_file = io.StringIO(bytes_file.decode("utf-8"))
        yield import_part(memory_file, "genbank")
