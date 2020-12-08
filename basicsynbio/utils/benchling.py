"""Module contains utilities to improve the interface with benchling."""


def only_label_feature(part):
    """Return part with only the label feature.

    Useful for getting rid of ApEinfo features following Benchling
    exports.

    Args:
        part: part subject to qualifier simplification.

    Returns:
        part: with simplified qualifiers
    """
    for feature in part.features:
        feature.qualifiers = {"label": feature.qualifiers["label"]}
    return part
