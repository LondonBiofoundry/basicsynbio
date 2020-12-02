from basicsynbio.decorators import add2docs
from basicsynbio.main import CommonArgDocs
from basicsynbio.main import BasicPart
from basicsynbio.cam import _seqrecord_hexdigest


class PartLinkerCollection(dict):
    def __str__(self):
        string = ""
        for item in self.items():
            string += f"{item[0]} id: {item[1].id}\n"
            string += f"{item[0]} name: {item[1].name}\n"
            string += f"{item[0]} description: {item[1].description}\n\n"
        return string


@add2docs(CommonArgDocs.PARTS_LINKERS_ARGS)
def make_collection(*parts_linkers, keys=None):
    """Generates a PartLinkerCollection object.
    Args:
        keys -- if None, uses id attribute, otherwise user supplies iterable of keys corresponding to each part/linker.
    """
    parts_linkers_withID = map(processId_part,parts_linkers)
    if not keys:
        collection = {part_linker.id: part_linker for part_linker in parts_linkers_withID}
    else:
        collection = {key: value for key, value in zip(keys, parts_linkers_withID)}
    return PartLinkerCollection(collection.items())

def processId_part(part):
    setattr(part,'id',_seqrecord_hexdigest(part))
    return part
