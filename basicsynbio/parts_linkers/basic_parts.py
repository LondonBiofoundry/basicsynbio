"""Module contains objects to access BASIC DNA assembly parts from within the pacakge."""

import basicsynbio as bsb
from importlib import resources
from .main import PartLinkerCollection, make_collection

with resources.open_text("basicsynbio.parts_linkers", f"BASIC_promoter_collection.gb") as gb_file:
    BPROMOTER_PARTS = make_collection(*bsb.import_parts(gb_file, "genbank"))

with resources.open_text("basicsynbio.parts_linkers", f"BASIC_CDS_collection.gb") as gb_file:
    BCDS_PARTS = make_collection(*bsb.import_parts(gb_file, "genbank"))

with resources.open_text("basicsynbio.parts_linkers", f"BASIC_SEVA_collection.gb") as gb_file:
    BSEVA_PARTS = list(bsb.import_parts(gb_file, "genbank"))
    BSEVA_PARTS = make_collection(
        *BSEVA_PARTS,
        keys=[part.id[-2:] for part in BSEVA_PARTS]
    )