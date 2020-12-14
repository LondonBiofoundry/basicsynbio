from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs, BasicPart, BasicLinker
from typing import Union, Iterable, Dict


class PartLinkerCollection(dict):
    """Class to hold Part Linker Collections.
    """
    def __str__(self):
        string = ""
        for item in self.items():
            string += f"{item[0]} id: {item[1].id}\n"
            string += f"{item[0]} name: {item[1].name}\n"
            string += f"{item[0]} description: {item[1].description}\n\n"
        return string


@add2docs(CommonArgDocs.PARTS_LINKERS_ARGS)
def make_collection(*parts_linkers: Union[BasicPart, BasicLinker], keys: Iterable[str] =None) -> Dict[str, Union[BasicPart, BasicLinker]]:
    """Generates a PartLinkerCollection object.

    Args:
        *parts_linkers (iterable of Parts/Linkers): iterable of BasicParts
            or BasicLinkers used to create the collection
        keys (optional): If None, uses id attribute, otherwise user supplies
            iterable of keys corresponding to each part/linker. Defaults to None

    Returns:
        Collection
            
    """
    if not keys:
        collection = {part_linker.id: part_linker for part_linker in parts_linkers}
    else:
        collection = {key: value for key, value in zip(keys, parts_linkers)}
    return PartLinkerCollection(collection.items())
