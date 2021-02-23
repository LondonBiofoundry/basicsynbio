from reportlab.platypus import (
    SimpleDocTemplate,
    ListFlowable,
    Paragraph,
    Spacer,
    Image,
    PageBreak,
    KeepTogether,
)
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.pdfgen import canvas
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.units import mm, cm
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT, TA_CENTER
from datetime import datetime


class PageNumCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        """Constructor"""
        canvas.Canvas.__init__(self, *args, **kwargs)
        self.pages = []

    def showPage(self):
        """
        On a page break, add information to the list
        """
        self.pages.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        """
        Add the page number to each page (page x of y)
        """
        page_count = len(self.pages)

        for page in self.pages:
            self.__dict__.update(page)
            self.draw_page_number(page_count)
            canvas.Canvas.showPage(self)

        canvas.Canvas.save(self)

    def draw_page_number(self, page_count):
        """
        Add the page number
        """
        page = "Page %s of %s" % (self._pageNumber, page_count)
        Title = f"Protocol for BasicBuild {datetime.now().strftime('%d-%m-%Y')}"
        self.setFont("Helvetica-Bold", 16)
        self.drawCentredString(A4[0] / 2.0, A4[1] - 2 * cm, Title)
        self.setFont("Helvetica", 9)
        self.drawCentredString(
            A4[0] / 2.0,
            1.1 * cm,
            "Page %d of %d" % (self._pageNumber, page_count),
        )


MATERIALS = [
    ["Item", "Order number (if applicable)"],
    ["PCR machine", ""],
    ["Water bath (42°C) for transformation", ""],
    ["Magnetic plate", "Ambion AM10050 (Thermo)"],
    ["Benchtop centrifuge", ""],
    ["Microplate centrifuge ", ""],
    ["Vortex", ""],
    ["96 well U-Bottom Falcon plate", "Falcon 351177 (Thermo)"],
    ["Microcentrifuge tubes", ""],
    ["PCR tubes", ""],
    ["Brooks Life Sciences PCR Foil Seals", "4ti-0550 "],
    ["Pipettes/tips 10 and 200 μl", ""],
    ["Agencourt AMPure XP SPRI paramagnetic beads", "A6388 (Beckman Coulter) "],
    ["ddH20", ""],
    ["70% EtOH", ""],
    ["Biolegio BASIC linkers", ""],
    ["NEB BsaI-HF v2 enzyme (R3733) 20 U/μl", "R3733 (NEB)"],
    ["Promega T4 ligase (M1801) 1-3U/μl", "M1801 (Promega)"],
    [
        "[optional] 10x Assembly Buffer: 0.2 M Tris:HCl (pH 8.0), 0.1 M MgCl2, 0.5 M KCl",
        "",
    ],
    [
        "Chemically competent cells (DH5alpha, 1x109 CFU/μg pUC19",
        "C2987I (NEB) or equivalent ",
    ],
    ["SOC media", ""],
    ["Petri dishes", ""],
    ["LB-Agar + antibiotic/s", ""],
]

INCUBATION_TIMES = [
    ["Temperature (°C)", "Time", ""],
    ["95", "2 min", ""],
    ["94 (-1 °C/cycle)", "40 seconds", "x70 cycles"],
    ["4", "Hold", ""],
]

TEMPERATURE_PROFILE_1 = [
    ["Temperature (°C)", "Time", ""],
    ["37", "2 min", "x20 cycles"],
    ["20", "1 min", ""],
    ["37", "5 min", ""],
]

ASSEMBLY_REACTION_3 = [
    ["Reagent", "Volume "],
    ["10x Assembly Buffer (or 10x NEB CutSmart)", "1 μl"],
    [
        "Each purified clip reaction required for the assembly",
        "1 μl for each",
    ],
    ["dH20 ", "Top up to 10 μl total volume"],
]

PCR_MACHINE_3 = [
    [
        "Temperature (°C)",
        "Time",
    ],
    [
        "50",
        "45 min",
    ],
    [
        "4",
        "Hold",
    ],
]

styles = getSampleStyleSheet()
styleN = styles["BodyText"]
styleN.alignment = TA_LEFT

style = TableStyle(
    [
        ("BACKGROUND", (0, 0), (-1, 0), colors.lightblue),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, -1), "Courier-Bold"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 5),
        ("LINEABOVE", (0, 0), (-1, -1), 1, colors.black),
        ("LINEBELOW", (0, 0), (-1, -1), 1, colors.black),
        ("LINEBEFORE", (0, 0), (-1, -1), 1, colors.black),
        ("LINEAFTER", (0, 0), (-1, -1), 1, colors.black),
    ]
)

style_no_split = TableStyle(
    [
        ("BACKGROUND", (0, 0), (-1, 0), colors.lightblue),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, -1), "Courier-Bold"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 5),
        ("LINEABOVE", (0, 0), (-1, -1), 1, colors.black),
        ("LINEBELOW", (0, 0), (-1, -1), 1, colors.black),
        ("LINEBEFORE", (0, 0), (-1, -1), 1, colors.black),
        ("LINEAFTER", (0, 0), (-1, -1), 1, colors.black),
        ("NOSPLIT", (0, 0), (-1, -1)),
    ]
)

style_joint_cell = TableStyle(
    [
        ("BACKGROUND", (0, 0), (-1, 0), colors.lightblue),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.whitesmoke),
        ("FONTNAME", (0, 0), (-1, -1), "Courier-Bold"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 5),
        ("LINEABOVE", (0, 0), (-1, -1), 1, colors.black),
        ("LINEBELOW", (0, 0), (-1, -1), 1, colors.black),
        ("LINEBELOW", (0, 0), (-1, -1), 0, colors.black),
        ("LINEBEFORE", (0, 0), (-1, -1), 1, colors.black),
        ("LINEAFTER", (0, 0), (-1, -1), 1, colors.black),
        ("SPAN", (2, 1), (2, 2)),
        ("NOSPLIT", (0, 0), (-1, -1)),
    ]
)

PROCESSED_MATERIALS = [list(map(lambda x: Paragraph(x, styleN), x)) for x in MATERIALS]
PROCESSED_INCUBATION_TIMES = [
    list(map(lambda x: Paragraph(x, styleN), x)) for x in INCUBATION_TIMES
]
PROCESSED_TEMPERATURE_PROFILE_1 = [
    list(map(lambda x: Paragraph(x, styleN), x)) for x in TEMPERATURE_PROFILE_1
]
PROCESSED_ASSEMBLY_REACTION_3 = [
    list(map(lambda x: Paragraph(x, styleN), x)) for x in ASSEMBLY_REACTION_3
]
PROCESSED_PCR_MACHINE_3 = [
    list(map(lambda x: Paragraph(x, styleN), x)) for x in PCR_MACHINE_3
]


def clips_data_from_pandas(clips_dataframe):
    headers = list(clips_dataframe.columns)
    del headers[2]
    data_for_clips_table = [headers]
    for row in clips_dataframe.to_numpy().tolist():
        del row[2]
        data_for_clips_table.append(row)
    return [
        list(map(lambda x: Paragraph(str(x), styleN), x)) for x in data_for_clips_table
    ]


def assembly_data_from_pandas(assembly_dataframe):
    data_for_assemblies_table = [list(assembly_dataframe.columns)]
    for row in assembly_dataframe.to_numpy().tolist():
        data_for_assemblies_table.append(row)
    return [
        list(map(lambda x: Paragraph(str(x), styleN), x))
        for x in data_for_assemblies_table
    ]
