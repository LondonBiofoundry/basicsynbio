import csv
from Bio.Seq import Seq
import json

biolegio_dict = {}
with open(r"csv_xlsx_files\Biolegio_linker_sequences.csv", newline="") as csv_file:
    reader = csv.reader(csv_file)
    for ind, line in enumerate(reader):
        if ind != 0:
            biolegio_dict[line[0]] = line[1].upper()

with open("biolegio.py", "w") as write_file:
    json.dump(biolegio_dict, write_file)
