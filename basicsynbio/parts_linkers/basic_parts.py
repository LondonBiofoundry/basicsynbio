"""Module contains objects to access BASIC DNA assembly parts from within the
pacakge."""

import basicsynbio as bsb
from importlib import resources
from .main import PartLinkerCollection, make_collection

with resources.open_text(
    "basicsynbio.parts_linkers", f"BASIC_promoter_collection.gb"
) as gb_file:
    imported_parts = list(bsb.import_parts(gb_file, "genbank"))
    BASIC_PROMOTER_PARTS = {
        "v0.1": make_collection(
            *imported_parts, keys=(part.id for part in imported_parts)
        )
    }

with resources.open_text(
    "basicsynbio.parts_linkers", f"BASIC_CDS_collection.gb"
) as gb_file:
    imported_parts = list(bsb.import_parts(gb_file, "genbank"))
    BASIC_CDS_PARTS = {
        "v0.1": make_collection(
            *imported_parts, keys=(part.id for part in imported_parts)
        )
    }

with resources.open_text(
    "basicsynbio.parts_linkers", f"BASIC_SEVA_collection.gb"
) as gb_file:
    BSEVA_PARTS = list(bsb.import_parts(gb_file, "genbank"))
    BASIC_SEVA_PARTS = {
        "v0.1": make_collection(
            *BSEVA_PARTS, keys=[part.id[-2:] for part in BSEVA_PARTS]
        )
    }
