"""Module contains objects to access BASIC DNA assembly parts from within the
pacakge."""

from basicsynbio.bsb_io import import_parts
from basicsynbio.main import (
    BasicPart,
    seqrec2part,
)
from importlib import resources
from basicsynbio.utils.primer3py import OLIGO_ANNOTATIONS
from .main import PartLinkerCollection, make_collection
from typing import (
    Generator,
    Iterable,
)


def parts_linkers_from_gb_file(
    file_name: str,
    package: str = "basicsynbio.parts_linkers",
    format: str = "genbank",
) -> Generator:
    """Import parts or linkers from a genbank file."""
    with resources.open_text(package, file_name) as gb_file:
        yield from import_parts(gb_file, format)


BASIC_CDS_PARTS_V01 = list(parts_linkers_from_gb_file("BASIC_CDS_collection_v01.gb"))
BASIC_CDS_PARTS_V02 = list(parts_linkers_from_gb_file("BASIC_CDS_collection_v02.gb"))
BASIC_CDS_PARTS = {
    "v0.1": make_collection(
        *BASIC_CDS_PARTS_V01, keys=(part.id for part in BASIC_CDS_PARTS_V01)
    ),
    "v0.2": make_collection(
        *BASIC_CDS_PARTS_V02, keys=(part.id for part in BASIC_CDS_PARTS_V02)
    ),
}

BASIC_PROMOTER_PARTS_V01 = list(
    parts_linkers_from_gb_file("BASIC_promoter_collection_v01.gb")
)
BASIC_PROMOTER_PARTS_V02 = list(
    parts_linkers_from_gb_file("BASIC_promoter_collection_v02.gb")
)
BASIC_PROMOTER_PARTS_V03 = list(
    parts_linkers_from_gb_file("BASIC_promoter_collection_v03.gb")
)
BASIC_PROMOTER_PARTS = {
    "v0.1": make_collection(
        *BASIC_PROMOTER_PARTS_V01, keys=(part.id for part in BASIC_PROMOTER_PARTS_V01)
    ),
    "v0.2": make_collection(
        *BASIC_PROMOTER_PARTS_V02, keys=(part.id for part in BASIC_PROMOTER_PARTS_V02)
    ),
    "v0.3": make_collection(
        *BASIC_PROMOTER_PARTS_V03, keys=(part.id for part in BASIC_PROMOTER_PARTS_V03)
    ),
}

BASIC_SEVA_PARTS_V01 = list(parts_linkers_from_gb_file("BASIC_SEVA_collection_v01.gb"))
BASIC_SEVA_PARTS_V10 = list(parts_linkers_from_gb_file("BASIC_SEVA_collection_v10.gb"))
BASIC_SEVA_PARTS_V20_DIFF = list(
    parts_linkers_from_gb_file("BASIC_SEVA_66-11_&_69-11.gb")
)
BASIC_SEVA_PARTS_V20 = BASIC_SEVA_PARTS_V10.copy()
BASIC_SEVA_PARTS_V20.pop(
    [part.name for part in BASIC_SEVA_PARTS_V20].index("BASIC_SEVA_66.10")
)
BASIC_SEVA_PARTS_V20.pop(
    [part.name for part in BASIC_SEVA_PARTS_V20].index("BASIC_SEVA_69.10")
)
BASIC_SEVA_PARTS_V20 += BASIC_SEVA_PARTS_V20_DIFF

BASIC_SEVA_PARTS_V10_KEYS = {
    "BASIC_SEVA_15a.10": "15a",
    "BASIC_SEVA_16.10": "16",
    "BASIC_SEVA_17.10": "17",
    "BASIC_SEVA_19.10": "19",
    "BASIC_SEVA_17_pKD46.10": "17_pKD46",
    "BASIC_SEVA_25a.10": "25a",
    "BASIC_SEVA_26.10": "26",
    "BASIC_SEVA_27.10": "27",
    "BASIC_SEVA_29.10": "29",
    "BASIC_SEVA_27_pKD46.10": "27_pKD46",
    "BASIC_SEVA_35a.10": "35a",
    "BASIC_SEVA_36.10": "36",
    "BASIC_SEVA_37.10": "37",
    "BASIC_SEVA_39.10": "39",
    "BASIC_SEVA_37_pKD46.10": "37_pKD46",
    "BASIC_SEVA_45a.10": "45a",
    "BASIC_SEVA_46.10": "46",
    "BASIC_SEVA_47.10": "47",
    "BASIC_SEVA_49.10": "49",
    "BASIC_SEVA_47_pKD46.10": "47_pKD46",
    "BASIC_SEVA_5a5a.10": "5a5a",
    "BASIC_SEVA_5a6.10": "5a6",
    "BASIC_SEVA_5a7.10": "5a7",
    "BASIC_SEVA_5a9.10": "5a9",
    "BASIC_SEVA_5a7_pKD46.10": "5a7_pKD46",
    "BASIC_SEVA_65a.10": "65a",
    "BASIC_SEVA_66.10": "66",
    "BASIC_SEVA_67.10": "67",
    "BASIC_SEVA_69.10": "69",
    "BASIC_SEVA_67_pKD46.10": "67_pKD46",
}
BASIC_SEVA_PARTS_V20_KEYS = BASIC_SEVA_PARTS_V10_KEYS.copy()
BASIC_SEVA_PARTS_V20_KEYS.pop("BASIC_SEVA_66.10")
BASIC_SEVA_PARTS_V20_KEYS.pop("BASIC_SEVA_69.10")
BASIC_SEVA_PARTS_V20_KEYS.update(
    {
        "BASIC_SEVA_66.11": "66",
        "BASIC_SEVA_69.11": "69",
    }
)
BASIC_SEVA_PARTS = {
    "v0.1": make_collection(
        *BASIC_SEVA_PARTS_V01, keys=(part.id[-2:] for part in BASIC_SEVA_PARTS_V01)
    ),
    "v1.0": make_collection(
        *BASIC_SEVA_PARTS_V10,
        keys=(BASIC_SEVA_PARTS_V10_KEYS[part.name] for part in BASIC_SEVA_PARTS_V10)
    ),
    "v2.0": make_collection(
        *BASIC_SEVA_PARTS_V20,
        keys=(BASIC_SEVA_PARTS_V20_KEYS[part.name] for part in BASIC_SEVA_PARTS_V20)
    ),
}
