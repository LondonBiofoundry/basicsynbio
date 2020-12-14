import basicsynbio as bsb
from typing import Union, Iterable, Iterator

myparts = bsb.import_parts('initial_BASIC_promoters.gb','genbank')
print(myparts)