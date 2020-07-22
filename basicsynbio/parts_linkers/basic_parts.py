"""Module contains objects to access BASIC DNA assembly parts from within the pacakge."""

import basicsynbio as bsb
from importlib import resources

BSEVA_DICT = {}
with resources.open_text("basicsynbio.parts_linkers", f"BASIC_SEVA_collection.gb") as gb_file:
    BSEVA_collection = bsb.import_parts(gb_file, "genbank")
    for part in BSEVA_collection:
        BSEVA_DICT[part.id[-2:]] = part
        
