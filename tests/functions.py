from basicsynbio.parts_linkers.main import PartLinkerCollection
from basicsynbio.cam import seqrecord_hexdigest
from typing import Union
import basicsynbio as bsb


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


def analyse_part_linker_collection(
    collection: PartLinkerCollection,
    collection_type: Union[bsb.BasicPart, bsb.BasicLinker],
    length: int,
):
    """Analyse a PartLinkerCollection instances to determine if it fufils basic requirements.

    Args:
        collection: PartLinkerCollection to analyse
        collection_type: collection contains either parts or linkers
        length: number of items in collection
    """
    assert len(collection) == length
    for part_linker in collection.values():
        assert part_linker.id == seqrecord_hexdigest(part_linker)
        assert type(part_linker) == collection_type
