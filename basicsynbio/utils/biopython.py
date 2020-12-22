"""Module contains utilities to improve the interface with Biopython."""

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from typing import Union, List


def feature_from_qualifier(
    seqrec: Union["BasicPart", SeqRecord], qualifier_key: str, qualifier_value: str
) -> SeqFeature:
    """Extract the feature from the seqrec that contains the corresponding
    qualifier key/qualifier value pair.

    Args:
        seqrec: sequence record subject to search.
        qualifier_key: search key within seqrec.features.
        qualifier_value: Search value within seqrec.features.

    Returns:
        feature: located feature with a key and value specified in Args.

    Raises:
        ValueError: If qualifier_key and qualifier_value cannot be found within
            provided seqrec.

    """
    for feature in seqrec.features:
        if qualifier_key in feature.qualifiers:
            if feature.qualifiers[qualifier_key] == qualifier_value:
                return feature
    raise ValueError(
        f"{seqrec.id} lacks a feature containing a {qualifier_key}/{qualifier_value} pair in qualifiers"
    )


def _easy_seqrec(
    str_seq: str,
    id,
    annotation_type: str = "misc_feature",
    start=0,
    end=None,
    **qualifiers: list,
) -> SeqRecord:
    """Return an annotated SeqRecord from a string and id.

    Args:
        str_seq: sequence of SeqRecord.
        id : Identifier for new part.
        annotation_type (optional): Equivalent to Bio.SeqFeature type e.g. CDS,
            Defaults to "misc_feature".
        start (optional): start of the annotation, Defaults to 0.
        end (optional): end of the annotation, if None defaults to len(str_seq),
            Defaults to 'None'.
        **qualifiers: equivalent to Bio.SeqFeature.qualifiers for annotation e.g. standard_name=["LMP"].

    Return:
        SeqRecord: An annotated SeqRecord.
    """
    if not end:
        end = len(str_seq)
    seqrec = SeqRecord(
        Seq(str_seq),
        id=id,
        features=[
            SeqFeature(
                type=annotation_type,
                location=FeatureLocation(start=start, end=end, strand=+1),
                qualifiers={item[0]: item[1] for item in qualifiers.items()},
            )
        ],
    )
    return seqrec


def all_feature_values(part: Union[SeqRecord, "BasicPart"]) -> List[str]:
    """Returns the values of all Bio.SeqFeature.qualifiers in part as a
    list.

    Args:
        part: part to be analysed

    Returns:
        list: the values of all Bio.SeqFeature.qualifiers
    """
    values = []
    for seqfeature in part.features:
        feature_values = list(seqfeature.qualifiers.values())
        values += [value for sublist in feature_values for value in sublist]
    return values
