from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


def feature_from_qualifier(seqrec, qualifier_key, qualifier_value):
    """Extract the feature from the seqrec that contains the corresponding qualifier key/qualifier value pair.

    """
    for feature in seqrec.features:
        if qualifier_key in feature.qualifiers:
            if feature.qualifiers[qualifier_key] == qualifier_value:
                return feature
    raise ValueError(f"{seqrec.id} lacks a feature containing a {qualifier_key}/{qualifier_value} pair in qualifiers")


def _easy_seqrec(str_seq: str, id, annotation_type="misc_feature", start=0, end=None,**qualifiers: list):
    """Return an annotated SeqRecord from a string and id.
        
    Args:
        str_seq -- sequence of SeqRecord
        id -- identifier for new part.
        annotation_type -- equivalent to Bio.SeqFeature type e.g. CDS
        start -- start of the annotation
        end -- end of the annotation, if None defaults to len(str_seq)
        **qualifiers -- equivalent to Bio.SeqFeature.qualifiers for annotation e.g. standard_name=["LMP"]
    
    """
    if not end:
        end = len(str_seq)
    seqrec = SeqRecord(Seq(str_seq), id=id, features=[SeqFeature(
        type=annotation_type,
        location=FeatureLocation(start=start, end=end, strand=+1),
        qualifiers={item[0]: item[1] for item in qualifiers.items()}
    )])
    return seqrec