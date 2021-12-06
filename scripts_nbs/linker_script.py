import csv
import pandas as pd
from pathlib import Path
import json

path_to_csv_xlsx = Path.cwd() / "csv_xlsx_files"
# dump dict for biolegio linkers
biolegio_dict = {}
with open(path_to_csv_xlsx / "Biolegio_linker_sequences.csv", newline="") as csv_file:
    reader = csv.reader(csv_file)
    for ind, line in enumerate(reader):
        if ind != 0:
            biolegio_dict[line[0]] = line[1].upper()
# with open("biolegio.py", "w") as write_file:
#     json.dump(biolegio_dict, write_file)
# dump dict for linker 96 plate layout
LINKER_PLATE_LAYOUT = pd.read_excel(
    path_to_csv_xlsx / "BASIC_linker_Biolegio_96well_update191128.xlsx",
    sheet_name="plate layout linkers",
    usecols="A:B",
)
LINKER_PLATE_LAYOUT = {
    row[1]["Well"]: row[1]["Linker"] for row in LINKER_PLATE_LAYOUT.iterrows()
}
with open("linker_plate_layout.txt", "w") as write_file:
    json.dump(LINKER_PLATE_LAYOUT, write_file, indent=4)
