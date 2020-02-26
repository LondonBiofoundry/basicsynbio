from Bio.SeqRecord import SeqRecord
from Bio import SeqUtils, SeqIO
from basicsynbio import basicsynbio_exceptions


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
        self.kwargs = kwargs

    def basic_slice(self):
        returned_seqrec = SeqRecord(seq=self.seq, id=self.id, **self.kwargs)
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
            raise basicsynbio_exceptions.PartException(f"{self.id} lacks {iseq_id}")
        return search_out[1]


def import_basic_part(handle, format, alphabet=None):
    """
    Arguments follow that of SeqIO.read
    """
    part = SeqIO.read(handle, format, alphabet)
    basic_part = BasicPart(part.seq, part.id)
    for key, value in part.__dict__.items():
        setattr(basic_part, key, value)
    return basic_part

