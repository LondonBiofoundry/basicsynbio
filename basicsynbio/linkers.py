from basicsynbio import BasicLinker
import csv
from Bio.Seq import Seq

biolegio_dict = {}
with open(r"csv_xlsx_files\Biolegio_linker_sequences.csv", newline='') as csv_file:
    reader = csv.reader(csv_file)
    for ind, line in enumerate(reader):
        if ind != 0:
            basic_linker = BasicLinker(Seq("GG" + line[1].upper()), line[0])
            biolegio_dict[line[0]] = basic_linker

if __name__ == "__main__":
    for key, value in biolegio_dict.items():
        print(f"{key}: {value}")