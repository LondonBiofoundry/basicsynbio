from basicsynbio.decorators import addargs2docs
from basicsynbio.cam import seqrecord_hexdigest
from basicsynbio.main import CommonArgDocs, BasicPart, BasicLinker
from typing import Union, Iterable, Dict


class PartLinkerCollection(dict):
    """Class to hold Part Linker Collections."""

    def __str__(self):
        string = ""
        for item in self.items():
            string += f"{'PartLinkerCollection key:':30} '{item[0]}'\n"
            string += f"{'id:':30} {item[1].id}\n"
            string += f"{'name:':30} {item[1].name}\n"
            string += f"{'description:':30} {item[1].description}\n\n"
        return string


@addargs2docs(CommonArgDocs.PARTS_LINKERS_ARGS)
def make_collection(
    *parts_linkers: Union[BasicPart, BasicLinker],
    keys: Iterable[str] = None,
    id_function: callable = None,
) -> Dict[str, Union[BasicPart, BasicLinker]]:
    """Generates a PartLinkerCollection object using parts_linkers.
    Args:
        parts_linkers:
        keys (optional): If None, uses id attribute, otherwise user supplies
            iterable of keys corresponding to each part/linker. Defaults to None
        id_function: function to define id of objects. If none uses set_part_linker_id function.

    Returns:
        Collection
    """
    parts_linkers_w_id = map(set_part_linker_id, parts_linkers)
    if not keys:
        collection = {
            part_linker.name: part_linker for part_linker in parts_linkers_w_id
        }
    else:
        collection = {key: value for key, value in zip(keys, parts_linkers_w_id)}
    return PartLinkerCollection(collection.items())


def set_part_linker_id(part_linker):
    """Sets the id attribute of a part_linker using the output of seqrecord_hexdigest."""
    part_linker.id = seqrecord_hexdigest(part_linker)
    return part_linker
