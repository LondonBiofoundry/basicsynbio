"""Module contains several small functions that support or utilise key cam
objects but each is not critical to the bsb workflow.
"""


from Bio.Restriction import BsaI
from Bio.SeqUtils import molecular_weight
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from basicsynbio.main import BasicAssembly, BasicPart, BasicLinker
from .main import BasicBuild
from dataclasses import dataclass
import hashlib
from typing import Generator, Tuple, Union


@dataclass
class AssemblyDigest:
    """class to return data calculated by `build_digest` function."""

    assembly: BasicAssembly
    restriction_enzyme: str
    sequences: Tuple[Seq, ...]
    product_lengths: Tuple[int, ...]


def seqrecord_hexdigest(seqrecord_obj: Union[SeqRecord, BasicPart, BasicLinker]) -> str:
    """Returns an MD5 hash of a Bio.SeqRecord.SeqRecord-like object, using
    relevant attributes.

    Note:
        Care should be taken when using the return value of this function as an identifier/hash for seqrecord_obj given these objects are mutable.

    Args:
        seqrecord_obj: Bio.SeqRecord.SeqRecord-like object, containing relevant
            attributes

    Returns:
        MD5 hash
    """
    seqrec_hash = hashlib.md5(str(seqrecord_obj.seq).encode("UTF-8"))
    bytes_objs = [
        getattr(seqrecord_obj, attribute).encode("UTF-8")
        for attribute in ["name", "description"]
    ]
    for element in bytes_objs:
        seqrec_hash.update(element)
    return seqrec_hash.hexdigest()


def new_part_resuspension(
    part: BasicPart, mass: float, double_stranded: bool = True
) -> float:
    """Returns the volume of resuspension buffer (µL) required for a 75 nM
    solution of part, equivalent to 75 fmol/µL.

    Args:
        part: BasicPart object.
        mass: mass of synthesised part (ng).
        double_stranded (optional): True (default) indicates part is dsDNA.
    """
    return (
        (mass * 10 ** -9)
        / molecular_weight(part.seq, double_stranded=double_stranded)
        * 1
        / (75e-9)
        * 10 ** 6
    )


def build_digest(
    build: BasicBuild,
    restriction_enzyme: None = BsaI,
) -> Generator[AssemblyDigest, None, None]:
    """Conduct a diagnostic digest on basic_assemblies in build.

    Function uses a single Bio.Restriction enzyme to digest all BasicAssembly objects
    in build. The resulting AssemblyDigest dataclass contains data on sequences and
    product lengths generated from digest.

    Args:
        build: BasicBuild object.
        restriction_enzyme: A single restriction enzyme from Bio.Restriction

    """
    for assembly in build.basic_assemblies:
        assembly_digest = restriction_enzyme.catalyse(
            assembly.return_part().seq, linear=False
        )
        yield AssemblyDigest(
            assembly,
            restriction_enzyme.__name__,
            assembly_digest,
            tuple(len(seq) for seq in assembly_digest),
        )
