"""basicsynbio - An open-source Python API to facilitate BASIC DNA Assembly workflows

When importing objects from modules, list each alphabetically.

"""

__author__ = "LondonBiofoundry"
__author_email__ = "hainesm6@gmail.com"
__description__ = "An open-source Python package to facilitate BASIC DNA Assembly workflows"
__python_version__ = ">=3.8"
__project__ = "basicsynbio"
__url__ = "https://github.com/LondonBiofoundry/basicsynbio"
__version__ = "develop"
__author__ = "Matthew Haines <hainesm6@gmail.com>"

__classifiers__=[
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Natural Language :: English",
    "Operating System :: OS Independent",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

__project_urls__={
    "Documentation": "https://londonbiofoundry.github.io/basicsynbio/index.html",
    "Source": "https://github.com/LondonBiofoundry/basicsynbio",
}

__requires__=[
    "biopython>=1.78",
    "pandas",
    "platemap",
    "primer3-py",
    "python-Levenshtein",
    "reportlab",
    "sbol2", 
]

from .main import BasicAssembly, BasicLinker, BasicPart, BasicUTRRBSLinker, seqrec2part
from .bsb_io import (
    export_sequences_to_file,
    import_part,
    import_parts,
    import_sbol_parts,
)
from .parts_linkers import (
    BASIC_CDS_PARTS,
    BASIC_BIOLEGIO_LINKERS,
    BASIC_PROMOTER_PARTS,
    BASIC_SEVA_PARTS,
)
from .cam import (
    BasicBuild,
    BuildDecoder,
    BuildEncoder,
    build_digest,
    new_part_resuspension,
    export_echo_assembly,
    pdf_instructions,
)
