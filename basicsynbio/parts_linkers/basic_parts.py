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
BASIC_SEVA_PARTS_V10_KEYS = {
    "f695f7983bc926fe22312ceea133ece8": "15a",
    "6568bcea77ad3ed95be9cf305b393522": "16",
    "50e0aa1b55ef02db0cf62a84c4ecb719": "17",
    "314fdd9dca836243b4181ccd287e2dc0": "19",
    "c2ce6224a5d3a0001db46ce89a840642": "17_pKD46",
    "4683b5e8bb9714875a28a0880c88f904": "25a",
    "b484d417d0b514190db028347fa66fba": "26",
    "3b433d3fd1e298bff6d5d066a372d067": "27",
    "b0dbd8bc2d60e1a7aad51780fcdbffd8": "29",
    "016b73059a8fa9268d8b35d2bdfa96ed": "27_pKD46",
    "01bb230a5ebe20dc3196b1f760d42d94": "35a",
    "1413513c978a56f243171c4f4cf6aaac": "36",
    "eceac71fa048bd052921bf134322c381": "37",
    "927f54a3afd7bab22d2c2cc3e428bf4a": "39",
    "b89dedf6aac9bef578d7fd8ac362c845": "37_pKD46",
    "56da3cbe26311835d21578bba476ebad": "45a",
    "4d797ce1665dcb6dc7e73376928424e7": "46",
    "b23db40e89aa7f017b096e4158e729bc": "47",
    "804fbb49400cf5027a668feb1ec60c06": "49",
    "a388b01c3624b66b0f4f70447ef4b42c": "47_pKD46",
    "4dc7813eba6298720b8e47f618aa2a4b": "5a5a",
    "da102cb0396b5d4633d32fa99b1b6cdb": "5a6",
    "d9f4ac855df0c3fcafcab2ae7e872036": "5a7",
    "de15e534a36c66f3ebe7c835782cffa7": "5a9",
    "1319bafa14d98407fe18627921dba1c3": "5a7_pKD46",
    "ba9809f0bba102f479d902c3fb82b5f7": "65a",
    "9874d4b94b5ebed8ea99f69534b4656e": "66",
    "db12df1f8d044fdd49a6237418780dd1": "67",
    "7cce40252f64427685493882230c8bf1": "69",
    "4fe6ecc619de0d790bdb44aeef93df39": "67_pKD46",
}
BASIC_SEVA_PARTS = {
    "v0.1": make_collection(
        *BASIC_SEVA_PARTS_V01, keys=(part.id[-2:] for part in BASIC_SEVA_PARTS_V01)
    ),
    "v1.0": make_collection(
        *BASIC_SEVA_PARTS_V10,
        keys=(BASIC_SEVA_PARTS_V10_KEYS[part.id] for part in BASIC_SEVA_PARTS_V10)
    ),
}
