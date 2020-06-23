from basicsynbio.utils import _easy_seqrec
from basicsynbio.decorators import add2docs
from Bio import SeqUtils, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils.CheckSum import seguid
import datetime

DATE = datetime.datetime.now()
DEFAULT_ANNOTATIONS = {
                        "source": "synthetic construct",
                        "organism": "synthetic construct",
                        "taxonomy": ["other sequences", "artificial sequences"],
                        "date": DATE.strftime("%d-") + DATE.strftime("%b").upper() + DATE.strftime("-%Y"),
                        "accessions": [],
                        "sequence_version": 1,
                        "topology": "circular"
                    }
IP_STR = "TCTGGTGGGTCTCTGTCC"
IS_STR = "GGCTCGGGAGACCTATCG"
IP_SEQREC = _easy_seqrec(IP_STR, "iP", note=["BASIC integrated prefix"])
IS_SEQREC = _easy_seqrec(IS_STR, "iS", note=["BASIC integrated suffix"])

class CommonArgDocs:
    ADD_I_SEQS = "add_i_seqs -- if True adds flanking BASIC iP and iS sequences. Note, letter_annotations attribute is lost."
    HANDLE = "handle -- handle to file."
    FORMAT = "format -- file format."
    ALPHABET = "alphabet -- Bio.Alphabet. Refer Bio.Alphabet documentation."
    SEQREC_KWARGS = "**kwargs -- assigns alternative Bio.SeqRecord attributes."        


class BasicPart(SeqRecord):
    """A DNA sequence that can be used in a BASIC DNA assembly.

    All sequences must contain intergated prefix and suffix sequences.
    """

    def __init__(self, seq, id, **kwargs):
        super().__init__(seq=seq, id=id, **kwargs)
        self._ip_loc = self._find_iseq(
            IP_STR, "iP sequence"
        )
        self._is_loc = self._find_iseq(
            IS_STR, "iS sequence"
        )

    def basic_slice(self):
        """Return the SeqRecord flanked by BASIC iP & iS sequences."""
        returned_seqrec = SeqRecord(seq=self.seq, id=self.id)
        for key in returned_seqrec.__dict__.keys():
            setattr(returned_seqrec, key, self.__dict__[key])
        if self._ip_loc < self._is_loc:
            return returned_seqrec[
                self._ip_loc + len(IP_STR):self._is_loc]
        elif self._ip_loc > self._is_loc:
            return returned_seqrec[self._ip_loc + len(IP_STR):] + returned_seqrec[:self._is_loc]
        else:
            raise ValueError("incorrect sequence used.")

    def _find_iseq(self, iseq_str, iseq_id="integrated sequence"):
        search_out = SeqUtils.nt_search(
            str(self.seq), iseq_str
        )
        if len(search_out) < 2:
            raise PartException(f"{self.id} lacks {iseq_id}")
        return search_out[1]


class BasicLinker(SeqRecord):
    def __init__(self, seq, id, **kwargs):
        super().__init__(seq=seq, id=id, **kwargs)
        self._linker_feature()

    def basic_slice(self):
        return self

    def _linker_feature(self):
        self.features.append(
            SeqFeature(
                type="misc_feature",
                location=FeatureLocation(2, len(self.seq), strand=+1),
                qualifiers={
                    "function": ["BASIC DNA assembly linker"],
                    "standard_name": [str(self.id)],
                    "note": [str(self.id)]
                }
            )
        )


class BasicAssembly():
    def __init__(self, *parts_linkers):
        """BasicAssembly class requires alternating BasicParts and BasicLinkers in any order.
        
        Args:
            *parts_linkers -- BasicPart or BasicLinker objects.
        """
        self.parts_linkers = parts_linkers

    @add2docs(
        12,
        CommonArgDocs.ALPHABET,
        CommonArgDocs.SEQREC_KWARGS
    )
    def return_part(self, id: str, alphabet=IUPAC.ambiguous_dna, **kwargs):
        """Return a new BASIC part from this assembly.

        Args:
            id -- identifier of the new part."""
        return seqrec2part(self.return_seqrec(alphabet=alphabet, id=id, **kwargs))
    
    @add2docs(
        12,
        CommonArgDocs.HANDLE,
        CommonArgDocs.FORMAT,
        CommonArgDocs.ALPHABET,
        CommonArgDocs.SEQREC_KWARGS
    )
    def return_file(self, handle, format="genbank", alphabet=IUPAC.ambiguous_dna, **kwargs):
        """Export BASIC assembly to "handle" using Bio.SeqIO.write().
        
        Args:"""
        seqrec = self.return_seqrec(alphabet, **kwargs)
        SeqIO.write(seqrec, handle, format)

    @add2docs(
        12,
        CommonArgDocs.ALPHABET,
        CommonArgDocs.SEQREC_KWARGS
    )
    def return_seqrec(self, alphabet=IUPAC.ambiguous_dna, **kwargs):
        """Returns a Bio.SeqRecord object of the assembled construct.
        
        Args:"""
        seqrec = SeqRecord(Seq(str()))
        for part_linker in self.parts_linkers:
            seqrec += part_linker.basic_slice()
        seqrec.id = seguid(seqrec.seq)
        seqrec.name = "BASIC_construct_" + seqrec.id
        seqrec.description = f"BASIC DNA Assembly of {[part_linker.id for part_linker in self.parts_linkers]}"
        seqrec.annotations = DEFAULT_ANNOTATIONS
        seqrec.seq.alphabet = alphabet
        if kwargs:
            for key, value in kwargs.items():
                setattr(seqrec, key, value)
        return seqrec

    @property
    def parts_linkers(self):
        return self._parts_linkers

    @parts_linkers.setter
    def parts_linkers(self, values):
        if not all(
                issubclass(type(value), BasicPart) or issubclass(type(value), BasicLinker) for value in values):
            raise TypeError(
                "Not all *parts_linkers are BasicParts or BasicLinkers."
            )
        for ind, value in enumerate(values):
            if ind != 0:
                if type(values[ind - 1]) == type(value):
                    raise AssemblyException(
                        f"Alternating BasicPart, BasicLinker instances required: {value.id} is preceeded by {type(values[ind - 1])} and is of type {type(value)}.")
        self._parts_linkers = values


class PartException(Exception):
    pass


class AssemblyException(Exception):
    pass


@add2docs(4, CommonArgDocs.ADD_I_SEQS)
def seqrec2part(seqrec, add_i_seqs=False):
    """Convert a Bio.SeqRecord to a BasicPart, relevant attributes are maintained.

    Args:"""
    if add_i_seqs:
        new_seqrec = IP_SEQREC + seqrec + IS_SEQREC
        part = BasicPart(new_seqrec.seq, seqrec.id, features=new_seqrec.features)
        for key, value in seqrec.__dict__.items():
            if key not in ("_seq", "_per_letter_annotations", "features"):
                setattr(part, key, value)
    else:
        part = BasicPart(seqrec.seq, seqrec.id)
        for key, value in seqrec.__dict__.items():
            setattr(part, key, value)
    return part


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
def export_to_file(parts, handle, format="genbank"):
    """Exports BasicPart/s using Bio.SeqIO.write().
    
    Args:
        parts -- BasicPart objects to export to file handle."""
    SeqIO.write(parts, handle, format)
