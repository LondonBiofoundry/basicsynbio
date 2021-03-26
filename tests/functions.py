def compare_seqrec_instances(seqrec1, seqrec2):
    """
    returns true if seqrec1 has equivalent seqrec2 attributes.
    Ignores seqrec.features as contains SeqFeature objects.

    """
    for key, value in seqrec2.__dict__.items():
        if key not in ["features", "basic_slice"]:
            if value != getattr(seqrec1, key):
                print(
                    f"seqrec2's '{key}' attribute does not match that obtained from seqrec1."
                )
                return False
    if len(seqrec1.features) != len(seqrec2.features):
        print("Length of features attribute is inconsistent between objects compared.")
        return False
    return True


def json_round_trip(class_instance, encoder, decoder):
    """Encodes object as serialised json and returns decoded object.

    Args:
        object -- class instance to be serialised.
        encoder -- json encoder Class.
        decoder -- json decoder function.
    """
    import json

    serialised_object = json.dumps(class_instance, cls=encoder)
    return json.loads(serialised_object, cls=decoder)
