"""Module contains utilities to improve the interface with benchling."""


def only_label_feature(part: "BasicPart") -> "BasicPart":
    """Return part with only the label feature.

    Useful for getting rid of ApEinfo features following Benchling
    exports.

    Args:
        part: object to simplify each feature in.

    Returns:
        part: with simplified qualifiers
    """
    for feature in part.features:
        feature.qualifiers = {"label": feature.qualifiers["label"]}
    return part
