"""Main module for basicsynbio."""

from basicsynbio.utils import _easy_seqrec
from basicsynbio.decorators import add2docs
from Bio import SeqUtils, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Restriction.Restriction import BsaI
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils.CheckSum import seguid
from collections import Counter
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
    PARTS_LINKERS_ARGS = "*parts_linkers -- BasicPart or BasicLinker objects."


class BasicPart(SeqRecord):
    """A DNA sequence that can be used in a BASIC DNA assembly.

    All sequences must contain intergated prefix and suffix sequences.
    """

    def __init__(self, seq, id, **kwargs):
        """Refer to help(BasicPart), id is required."""
        super().__init__(seq=seq, id=id, **kwargs)
        self._ip_loc = self._find_iseq(
            IP_STR, "iP sequence"
        )
        self._is_loc = self._find_iseq(
            IS_STR, "iS sequence"
        )
        self._check_bsai()

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
        elif len(search_out) > 2:
            raise PartException(f"{self.id} contains multiple {iseq_id}")
        return search_out[1]

    def _check_bsai(self):
        """Checks if sliced basic part contains a BsaI site."""
        if len(BsaI.search(self.seq)) > 2:
            raise PartException(
                f"{self.id} contains more than two BsaI sites.")

    def __eq__(self, other):
        if not isinstance(other, BasicPart):
            raise TypeError(f"{other} is not a BasicPart instance.")
        return (
            self.id == other.id and
            str(self.seq) == str(other.seq)
        )


class BasicLinker(SeqRecord):
    def __init__(self, seq, id, prefix_id=None, suffix_id=None, **kwargs):
        """Constructor for BasicLinker.

        Args:
            seq -- Refer to BioPython.SeqRecord.SeqRecord documentation.
            id -- Refer to BioPython.SeqRecord.SeqRecord documentation
            prefix_id -- ID for prefix linker half.
            suffix_id -- ID for suffix linker half.
        """
        super().__init__(seq=seq, id=id, **kwargs)
        self.prefix_id = self._assign_linker_half_id("prefix", prefix_id)
        self.suffix_id = self._assign_linker_half_id("suffix", suffix_id)
        self._linker_feature()

    def basic_slice(self):
        return self

    def _linker_feature(self):
        self.features.append(
            SeqFeature(
                type="misc_feature",
                location=FeatureLocation(2, len(self.seq), strand=+1),
                qualifiers={
                    "label": [str(self.id)]
                }
            )
        )

    def _assign_linker_half_id(self, linker_half, id):
        if not id and linker_half == "prefix":
            return f"{self.id}-P"
        elif not id and linker_half == "suffix":
            return f"{self.id}-S"
        return id

    def __eq__(self, other):
        if not isinstance(other, BasicLinker):
            raise TypeError(f"{other} is not a BasicLinker instance.")
        return (
            self.id == other.id and
            str(self.seq) == str(other.seq)
        )


class BasicUTRRBSLinker(BasicLinker):
    """Sub-class for UTR-RBS linkers."""

    def __init__(self, seq, id, prefix_id=None, suffix_id=None, **kwargs):
        super().__init__(seq, id, prefix_id, suffix_id, **kwargs)
        self.prefix_id = super()._assign_linker_half_id("prefix", prefix_id)
        self.suffix_id = f"UTR{self.id[3]}-S"


class BasicAssembly():
    @add2docs(
        12,
        CommonArgDocs.PARTS_LINKERS_ARGS
    )
    def __init__(self, id: str, *parts_linkers):
        """BasicAssembly class requires alternating BasicParts and BasicLinkers in any order.

        Args:
            id -- Identifier for BasicAssemby object. Should be unique amongst BasicAssembly instances in a BasicBuild."""
        self.id = id
        self.parts_linkers = parts_linkers
        self.clip_reactions = self.return_clip_reactions()

    @add2docs(
        12,
        CommonArgDocs.ALPHABET,
        CommonArgDocs.SEQREC_KWARGS
    )
    def return_part(self, alphabet=IUPAC.ambiguous_dna, **kwargs):
        """Return a new BASIC part from this assembly.

        Args:"""
        return seqrec2part(self.return_seqrec(alphabet=alphabet, **kwargs))

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
        seqrec.id = self.id
        seqrec.name = "BASIC_construct_" + self.id
        seqrec.description = f"BASIC DNA Assembly of {[part_linker.id for part_linker in self.parts_linkers]}"
        seqrec.annotations = DEFAULT_ANNOTATIONS
        seqrec.seq.alphabet = alphabet
        if kwargs:
            for key, value in kwargs.items():
                setattr(seqrec, key, value)
        return seqrec

    def return_clip_reactions(self):
        """returns clip reactions required for assembly of self."""
        clip_reactions = []
        for ind, part_linker in enumerate(self.parts_linkers):
            if issubclass(type(part_linker), BasicLinker):
                pass
            else:
                if ind == len(self.parts_linkers) - 1:
                    suffix = self.parts_linkers[0]
                else:
                    suffix = self.parts_linkers[ind+1]
                clip_reactions.append(
                    ClipReaction(
                        prefix=self.parts_linkers[ind-1],
                        part=part_linker,
                        suffix=suffix
                    )
                )
        self._check_clip_reactions(clip_reactions)
        return tuple(clip_reactions)

    def _check_clip_reactions(self, clip_reactions):
        """Checks clip reactions are compatible e.g. same half linker not used multiple times."""

        def _check_linker_halves(linker_halves):
            """Check linker_havles are compatible. Note UTR linker-halves must be compatible."""
            if len(linker_halves) > len(set(linker_halves)):
                linker_half_counts = Counter(linker_halves)
                raise AssemblyException(
                    f"BasicAssembly initiated with {linker_half_counts.most_common(1)[0][0]} used {linker_half_counts.most_common(1)[0][1]} times.")

        prefix_linkers = [clip_reaction.linker_half_ids()[0]
                          for clip_reaction in clip_reactions]
        _check_linker_halves(prefix_linkers)
        suffix_linkers = [clip_reaction.linker_half_ids()[1]
                          for clip_reaction in clip_reactions]
        _check_linker_halves(suffix_linkers)

    @property
    def parts_linkers(self):
        return self._parts_linkers

    @parts_linkers.setter
    def parts_linkers(self, values):
        if not all(
                issubclass(type(value), BasicPart) or issubclass(type(value), BasicLinker) for value in values):
            raise TypeError(
                "Not all *parts_linkers are BasicParts or BasicLinkers instances."
            )
        for ind, value in enumerate(values):
            if ind != 0:
                if type(values[ind - 1]) == type(value):
                    raise AssemblyException(
                        f"Alternating BasicPart, BasicLinker instances required: {value.id} is preceeded by {type(values[ind - 1])} and is of type {type(value)}.")
        self._parts_linkers = values


class ClipReaction():
    """Class for describing clip reactions. Note ClipReaction is hashable."""

    def __init__(self, prefix, part, suffix):
        """Initiaties ClipReaction.

        Args:
            prefix -- BasicLinker instance used as prefix in clip reaction.
            part -- BasicPart instance.
            suffix -- BasicLinker instance used as suffix in clip reaction.

        """
        self._prefix = prefix
        self._suffix = suffix
        self._part = part

    def linker_half_ids(self):
        """Returns the ids for prefix and suffix linkers in the form (prefix_id, suffix_id)."""
        return self._prefix.prefix_id, self._suffix.suffix_id

    def clip_items(self):
        """Return (prefix, part, suffix)."""
        return self._prefix, self._part, self._suffix

    def __hash__(self):
        return hash((self._prefix.id, self._part.id, self._suffix.id))

    def __eq__(self, other):
        if not isinstance(other, ClipReaction):
            raise TypeError(f"{other} is not a ClipReaction instance.")
        return (
            self._prefix == other._prefix and
            self._part == other._part and
            self._suffix == other._suffix
        )


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
        part = BasicPart(new_seqrec.seq, seqrec.id,
                         features=new_seqrec.features)
        for key, value in seqrec.__dict__.items():
            if key not in ("_seq", "_per_letter_annotations", "features"):
                setattr(part, key, value)
    else:
        part = BasicPart(seqrec.seq, seqrec.id)
        for key, value in seqrec.__dict__.items():
            setattr(part, key, value)
    return part
