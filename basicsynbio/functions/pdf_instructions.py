from basicsynbio.cam import BasicBuild
import zipfile
import os
import pandas as pd
from pathlib import Path
from datetime import datetime
from reportlab.platypus import (
    SimpleDocTemplate,
    ListFlowable,
    Paragraph,
    Spacer,
    Image,
    PageBreak,
    KeepTogether,
    Table,
)
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from basicsynbio.utils import (
    PROCESSED_MATERIALS,
    PROCESSED_INCUBATION_TIMES,
    PROCESSED_TEMPERATURE_PROFILE_1,
    PROCESSED_ASSEMBLY_REACTION_3,
    PROCESSED_PCR_MACHINE_3,
    PageNumCanvas,
    style,
    style_joint_cell,
    style_no_split,
    styleN,
    styles,
    clips_data_from_pandas,
    assembly_data_from_pandas,
)


def pdf_instructions(basic_build: BasicBuild, assemblies_per_clip: int = 28):
    """Writes information about each clip_data and assembly to
    two dependent CSV files in the same folder the command
    is executed.

    Args:
        path (optional): path to zipped folder of csv files. If none defaults to
            working directory with a time stamped name, output csvs is created.

    Returns:
        str: filepath of created pdf
    """

    import pandas as pd
    import math

    def calculate_clip_num(assemblies: list):
        return math.ceil(len(assemblies) / assemblies_per_clip)

    COMPONENTS = pd.DataFrame(
        {
            "Component": [
                "Promega T4 DNA Ligase 10x Buffer",
                "Water",
                "NEB BsaI-HFv2",
                "Promega T4 DNA Ligase",
            ],
            "Volume per clip (µL)": [3, 15.5, 1, 0.5],
        }
    )

    total_clips = sum(
        [
            calculate_clip_num(assemblies)
            for _, assemblies in basic_build.clips_data.items()
        ]
    )
    dead_clips = math.ceil(total_clips / 20)
    array = COMPONENTS["Volume per clip (µL)"] * (total_clips + dead_clips)

    MASTER_MIX_BASIC_REACTION = [
        [
            "Component",
            "Volume per clip (µL)",
        ],
        [
            "Promega T4 DNA Ligase 10x Buffer",
            str(array[0]),
        ],
        [
            "Water",
            str(array[1]),
        ],
        [
            "NEB BsaI-HFv2",
            str(array[2]),
        ],
        [
            "Promega T4 DNA Ligase",
            str(array[3]),
        ],
    ]

    PROCESSED_MASTER_MIX_BASIC_REACTION = [
        list(map(lambda x: Paragraph(x, styleN), x)) for x in MASTER_MIX_BASIC_REACTION
    ]

    zip_path = basic_build.export_csvs()
    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        try:
            zip_ref.extractall("PDF_CSVS")
        finally:
            zip_ref.close()
            os.remove(zip_path)
    assemeblies_data = pd.read_csv(Path.cwd() / "PDF_CSVS" / "assemblies.csv")
    clips_data = pd.read_csv(Path.cwd() / "PDF_CSVS" / "clips.csv")
    os.remove(Path.cwd() / "PDF_CSVS" / "assemblies.csv")
    os.remove(Path.cwd() / "PDF_CSVS" / "clips.csv")
    os.rmdir(Path.cwd() / "PDF_CSVS")
    CLIPS_DATA = clips_data_from_pandas(clips_data)
    ASSEMBLIES_DATA = assembly_data_from_pandas(assemeblies_data)

    pdf_filename = f"pdf_{datetime.now().strftime('%d-%m-%Y_%H.%M.%S')}.pdf"
    pdf = SimpleDocTemplate(pdf_filename, pagesizes=A4)

    elems = []
    elems.append(Image("basicsynbio/static/introimg.png", 12 * cm, 4 * cm))
    elems.append(Paragraph("Materials", styles["Heading1"]))
    elems.append(Table(PROCESSED_MATERIALS, colWidths=[10.5 * cm, 5 * cm], style=style))
    elems.append(PageBreak())
    elems.append(Paragraph("Method", styles["Heading1"]))
    elems.append(
        Paragraph("Preperation of BASIC linkers and parts", styles["Heading2"])
    )
    elems.append(
        Paragraph(
            "BASIC linkers can be ordered from Biolegio (info@biolegio.com) and will be delivered in a lyophilized format along with linker annealing buffer.",
            styles["BodyText"],
        )
    )
    elems.append(
        ListFlowable(
            [
                Paragraph(
                    "Spin down plates containing lyophilized linkers. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Add 150 μl of linker annealing buffer to each well, seal the plate with a PCR foil and incubate for 1 hour at room temperate. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Vortex the plate and collect liquid via centrifugation. ",
                    styles["BodyText"],
                ),
                KeepTogether(
                    [
                        Paragraph(
                            "Conduct the following incubation in a thermocycler: ",
                            styles["BodyText"],
                        ),
                        Spacer(1, 0.4 * cm),
                        Table(
                            PROCESSED_INCUBATION_TIMES,
                            colWidths=[4 * cm, 4 * cm, 4 * cm],
                            style=style_no_split,
                        ),
                    ]
                ),
                Paragraph(
                    "Collect the liquid in wells via centrifugation.",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Linkers are ready to use or can be stored at -20°C until required.",
                    styles["BodyText"],
                ),
            ]
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    ## Basic Reaction 1
    elems.append(
        Paragraph(
            "Clips Reaction",
            styles["Heading2"],
        )
    )
    elems.append(
        Paragraph(
            "Below contains a table with each clip needed for the assemblies of the build.",
            styles["BodyText"],
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    elems.append(
        Table(
            CLIPS_DATA,
            colWidths=[
                2 * cm,
                1.9 * cm,
                3 * cm,
                3 * cm,
                2 * cm,
                2 * cm,
                2.3 * cm,
                2 * cm,
            ],
            style=style,
        )
    )
    elems.append(
        Paragraph(
            "Prepare a Master mix for clip reactions, the below table provides the required components for the master mix with sufficient quantities for all clip reactions.",
            styles["BodyText"],
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    elems.append(
        Table(
            PROCESSED_MASTER_MIX_BASIC_REACTION,
            colWidths=[7 * cm, 7 * cm],
            style=style_no_split,
        )
    )
    elems.append(
        Paragraph(
            "For each Clip reaction, setup 1 PCR tube with 30 μl total volume: ",
            styles["BodyText"],
        )
    )
    elems.append(
        Paragraph(
            "Dispense 20 uL master mix, 1 μl of each prefix and suffix Linker, 1 μl of part (or more depending on concentration) into a PCR tube and make up to 30 μl with water.",
            styles["BodyText"],
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    elems.append(
        Paragraph(
            "We recommend the above up each clip can be used in up to 28 assemblies.",
            styles["BodyText"],
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    elems.append(
        Paragraph(
            "After mixing, tubes are placed in a PCR machine running the following programme: ",
            styles["BodyText"],
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    elems.append(
        Table(
            PROCESSED_TEMPERATURE_PROFILE_1,
            colWidths=[4 * cm, 4 * cm, 4 * cm],
            style=style_joint_cell,
        )
    )
    ## Magbead Purification 2
    elems.append(Paragraph("Magbead purification", styles["Heading2"]))
    elems.append(
        Paragraph(
            "Prepare fresh 70% EtOH (0.5 ml per BASIC reaction) and bring magnetic beads (AmpureXP or Ampliclean) stored at 4°C back into homogeneous mix by shaking thoroughly. ",
            styles["BodyText"],
        )
    )
    elems.append(
        ListFlowable(
            [
                Paragraph(
                    "Add 54 μl of magnetic beads into 96 well Falcon plate (one well per BASIC reaction) and add the 30 μl BASIC linker ligation from the PCR machine step, mix by pipetting 10 times. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Wait 5 min to allow DNA binding to magbeads. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Place Falcon plate on magnetic stand and wait for rings to form and solution to clear.",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Aspirate most of the solution with a 200 µL pipette set to 80 µL.",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Add 190 μl 70% EtOH to each well and wait 30 s. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Remove solution from each well (pipette set to 200 μl volume) ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Add 190 μl 70% EtOH to each well and wait 30 s. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Remove solution from each well (pipette set to 200 μl volume)",
                    styles["BodyText"],
                ),
                Paragraph("Leave the plate to dry for 1-2 min.", styles["BodyText"]),
                Paragraph(
                    "Remove Falcon plate from magnet and resuspend magbeads in 32 μl dH20. ",
                    styles["BodyText"],
                ),
                Paragraph("Wait 1 min for DNA to elute. ", styles["BodyText"]),
                Paragraph(
                    "Place Falcon plate back on magnetic stand and allow ring to form and solution to clear. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Transfer 30 μl of purified clip reaction into a clean microcentrifuge tube or well and store at at -20°C for up to 1 month.",
                    styles["BodyText"],
                ),
            ]
        )
    )
    # Assembly Reaction 3
    elems.append(
        Paragraph(
            "Assembly reaction",
            styles["Heading2"],
        )
    )
    elems.append(
        Paragraph(
            "Below contains a table with each BASIC assembly within your build",
            styles["BodyText"],
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    elems.append(
        Table(ASSEMBLIES_DATA, colWidths=[5.2 * cm, 5.2 * cm, 5.2 * cm], style=style)
    )
    elems.append(
        Paragraph(
            "For each BASIC assembly, combine the required purified clip reactions in a final volume of 10 μl in 1x Assembly or NEB CutSmart buffer. Below gives an example for a 3-part assembly:",
            styles["BodyText"],
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    elems.append(
        Table(
            PROCESSED_ASSEMBLY_REACTION_3,
            colWidths=[7 * cm, 4 * cm],
            style=style_no_split,
        )
    )
    elems.append(
        Paragraph(
            "Run assembly reaction in PCR machine with following programme",
            styles["BodyText"],
        )
    )
    elems.append(Spacer(1, 0.4 * cm))
    elems.append(
        Table(
            PROCESSED_PCR_MACHINE_3,
            colWidths=[7 * cm, 4 * cm],
            style=style_no_split,
        )
    )
    # Assembly Reaction 4
    elems.append(
        Paragraph(
            "Transformation",
            styles["Heading2"],
        )
    )
    elems.append(
        Paragraph(
            "Use 50 μl of chemically competent cells DH5alpha with high transformation efficiency  (109 CFU/μg pUC19, for instance NEB C2987I) to transform 5 μl of each BASIC assembly: ",
            styles["BodyText"],
        )
    )
    elems.append(
        ListFlowable(
            [
                Paragraph(
                    "Chemically competent cells are stored at -80°C. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Thaw competent cells on ice (takes 5-10 min); 50 μl per BASIC assembly to be transformed. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Cool 5 μl of BASIC DNA assembly in 1.5 ml microcentrifuge tube on ice. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Add 50 μl of competent cells to each precooled 5 μl BASIC reaction. ",
                    styles["BodyText"],
                ),
                Paragraph("Incubate on ice for 20 min. ", styles["BodyText"]),
                Paragraph(
                    "Apply heat shock in 42°C water bath for 45s and place back on ice for 2 min. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Add 200 μl SOC medium to each tube and incubate shaking at 37°C for 1h recovery. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Spot or plate cells on agar plates with appropriate antibiotics. Depending on number of parts assembled and transformation efficiency 2-250 μl might be spotted or plated. ",
                    styles["BodyText"],
                ),
                Paragraph(
                    "Incubate agar plates at 37°C overnight, next day pick colony for assay or miniprep. ",
                    styles["BodyText"],
                ),
            ]
        )
    )
    pdf.build(elems, canvasmaker=PageNumCanvas)

    return Path.cwd() / pdf_filename
