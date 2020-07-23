from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs


class PartLinkerCollection(dict):
    def __str__(self):
        string = ""
        for item in self.items():
            string += f"{item[0]} id: {item[1].id}\n"
            string += f"{item[0]} name: {item[1].name}\n"
            string += f"{item[0]} description: {item[1].description}\n\n"
        return string


@add2docs(
    8,
    CommonArgDocs.PARTS_LINKERS_ARGS
)
def make_collection(*parts_linkers, key_mapping=None):
    """Generates a PartLinkerCollection object.
    
    Args:
        key_mapping -- if None, uses id attribute, otherwise user supplies iterable of keys corresponding to each part/linker."""
    if not key_mapping:
        collection = {part_linker.id: part_linker for part_linker in parts_linkers}
    else:
        collection = dict()
        for ind, key in enumerate(key_mapping):
            collection[key] = parts_linkers[ind]
    return PartLinkerCollection(collection.items())