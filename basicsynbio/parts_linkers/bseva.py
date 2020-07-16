import basicsynbio as bsb
from importlib import resources

BSEVA_DICT = {}
for ind in range(12, 69):
    # skip SEVA genbank files with certain last digits
    if str(ind)[-1] in [str(num) for num in (0, 1, 6)]:
        pass
    else:
        with resources.open_text("basicsynbio.parts_linkers", f"BASIC_SEVA_{ind}.gb") as gb_file:
            BSEVA_DICT[str(ind)] = bsb.import_part(gb_file, "genbank")
