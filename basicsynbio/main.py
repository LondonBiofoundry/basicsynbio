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

class BasicPart(SeqRecord):
    IP_STR = "TCTGGTGGGTCTCTGTCC"
    IS_STR = "GGCTCGGGAGACCTATCG"

    def __init__(self, seq, id, **kwargs):
        super().__init__(seq=seq, id=id, **kwargs)
        self._ip_loc = self._find_iseq(
            self.IP_STR, "iP sequence"
        )
        self._is_loc = self._find_iseq(
            self.IS_STR, "iS sequence"
        )

    def basic_slice(self):
        returned_seqrec = SeqRecord(seq=self.seq, id=self.id)
        for key in returned_seqrec.__dict__.keys():
            setattr(returned_seqrec, key, self.__dict__[key])
        if self._ip_loc < self._is_loc:
            return returned_seqrec[
                self._ip_loc + len(self.IP_STR):self._is_loc]
        elif self._ip_loc > self._is_loc:
            return returned_seqrec[self._ip_loc + len(self.IP_STR):] + returned_seqrec[:self._is_loc]
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
                    "standard_name": [f"{self.id}"]
                }
            )
        )


class BasicAssembly():
    def __init__(self, *parts_linkers):
        """
        BasicAssembly class requires alternating BasicPart and BasicLinkers

        """
        self.parts_linkers = parts_linkers

    def _return_seqrec(self, **kwargs):
        seqrec = SeqRecord(Seq(str()))
        for part_linker in self.parts_linkers:
            seqrec += part_linker.basic_slice()
        seqrec.id = seguid(seqrec.seq)
        seqrec.name = "BASIC_construct_" + seqrec.id
        seqrec.description = f"BASIC DNA Assembly of {[part_linker.id for part_linker in self.parts_linkers]}"
        seqrec.annotations = DEFAULT_ANNOTATIONS
        if kwargs:
            for key, value in kwargs.items():
                setattr(seqrec, key, value)
        return seqrec

    def return_file(self, handle, format="genbank", alphabet=IUPAC.ambiguous_dna, **kwargs):
        """
        export BASIC assembly to "handle" using Bio.SeqIO.write().
        **kwargs assigns alternative Bio.SeqRecord attributes.

        """
        seqrec = self._return_seqrec(**kwargs)
        seqrec.seq.alphabet = alphabet
        SeqIO.write(seqrec, handle, format)

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


def import_part(handle, format):
    """
    returns a BASIC part object using Bio.SeqIO.read()

    """
    seqrec = SeqIO.read(handle, format)
    return _seqrec2part(seqrec)
    

def import_parts(handle, format):
    """
    returns BASIC parts in a single file using Bio.SeqIO.parse()

    """
    parts = SeqIO.parse(handle, format)
    return [_seqrec2part(part) for part in parts]


def _seqrec2part(seqrec):
    part = BasicPart(seqrec.seq, seqrec.id)
    for key, value in seqrec.__dict__.items():
        setattr(part, key, value)
    return part