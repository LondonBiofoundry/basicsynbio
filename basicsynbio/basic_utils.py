def feature_from_qualifier(seqrec, qualifier_key, qualifier_value):
    """
    Extract the feature from the seqrec that contains the corresponding
    qualifier key/qualifier value pair.

    """
    for feature in seqrec.features:
        if qualifier_key in feature.qualifiers:
            if feature.qualifiers[qualifier_key] == qualifier_value:
                return feature
    raise ValueError(f"{seqrec.id} lacks a feature containing a {qualifier_key}/{qualifier_value} pair in qualifiers")