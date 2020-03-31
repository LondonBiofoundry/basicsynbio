import basicsynbio as bsb
from basicsynbio.main import DEFAULT_ANNOTATIONS
from Bio import SeqIO, Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import nt_search
from dataclasses import dataclass
from pathlib import Path

MODULE_DIR = (Path.cwd() / "scripts_nbs" / "dnabot_genbank_files")

def main():
    generate_vectors()
    generate_parts()


def generate_vectors():
    paths = [
        "basic_seva_18_ampr-puc-1.gb",
        "basic_seva_37_cmr-p15a-1.gb"
    ]
    paths = [(MODULE_DIR / path) for path in paths]
    vectors = (bsb.import_part(handle, "genbank") for handle in paths)
    linkers_features = (bsb.biolegio_dict["LMP"], bsb.biolegio_dict["LMS"])
    for vector in vectors:
        vector.id = f"{vector.id}1"
        vector.description="Vector for BASIC DNA assembly containing \
            mScarlet counter selection marker"
        vector.name = vector.id[:len(vector.id)-2]
        vector.annotations = DEFAULT_ANNOTATIONS
        for linker in linkers_features:
            add_feature_w_seq(
                vector,
                str(linker.features[0].extract(linker.seq)),
                function=["BASIC DNA assembly linker"],
                standard_name=[str(linker.id)],
                note=[str(linker.id)]
            )
        SeqIO.write(vector, (MODULE_DIR / f"{vector.id}.gb"), "genbank")


@dataclass
class QuickBasicPart:
    seq: str
    abbreviation: str
    part_type: str


def generate_parts():
    def make_gb_from_quick_part(quick_basic_part):
        if quick_basic_part.part_type == "promoter":
            name = f"J23{quick_basic_part.abbreviation}"
            handle = f"BASIC_L3S2P21_J23{quick_basic_part.abbreviation}_RiboJ.1.gb"
            annotation_type = "regulatory"
            qualifiers = {
                "note": ["promoter"],
                "standard_name": [name]
            }
        elif quick_basic_part.part_type == "orf":
            handle = f"BASIC_{quick_basic_part.abbreviation}_ORF.1.gb"
            name = quick_basic_part.abbreviation
            annotation_type = "CDS"
            qualifiers = {
                "gene": [quick_basic_part.abbreviation],
                "note": ["fluorescent reporter protein"],
                "translation": [str(Seq.Seq(quick_basic_part.seq).translate()[:-1])]
            }
        else:
            raise ValueError(f"{quick_basic_part.abbreviation} part_type not compatible.")
        
        basic_part_creator = bsb.BasicPartCreator(
            str_seq=quick_basic_part.seq,
            id=name,
            annotation_type=annotation_type,
            **qualifiers
        )
        assembly = bsb.BasicAssembly(
            bsb.biolegio_dict["LMS"],
            bsb.import_part((MODULE_DIR / "BASIC_SEVA_18_AmpR-pUC.1.gb"), "genbank"),
            bsb.biolegio_dict["LMP"],
            basic_part_creator.create_part()
        )
        assembly.return_file(
            (MODULE_DIR / handle),
            id=handle[:len(handle)-3],
            name=handle[:len(handle)-5],
            annotations=DEFAULT_ANNOTATIONS,
            description= f"{name} BASIC part stored in AmpR pUC vector"
        )

    parts = []
    parts.append(QuickBasicPart(
        "CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACGGCTAGCTCAGTCCTAGGTACTATGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", "105", "promoter"
    ))
    parts.append(QuickBasicPart(
        "CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACGGCTAGCTCAGTCCTAGGTATAGTGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", "106", "promoter"
    ))
    parts.append(QuickBasicPart(
        "CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTTACAGCTAGCTCAGTCCTAGGTATTATGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", "101", "promoter"
    ))
    parts.append(QuickBasicPart(
        "CTCGGTACCAAATTCCAGAAAAGAGGCCTCCCGAAAGGGGGGCCTTTTTTCGTTTTGGTCCGTGCCTACTCTGGAAAATCTTTGACAGCTAGCTCAGTCCTAGGTATTGTGCTAGCAGCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGCCTCTACAAATAATTTTGTTTAA", "104", "promoter"
    ))
    parts.append(QuickBasicPart(
        "ATGCGTAAAGGCGAAGAACTGTTCACGGGCGTAGTTCCGATTCTGGTCGAGCTGGACGGCGATGTGAACGGTCATAAGTTTAGCGTTCGCGGTGAAGGTGAGGGCGACGCGACCAACGGCAAACTGACCCTGAAGTTCATCTGCACCACCGGTAAACTGCCGGTGCCTTGGCCGACCTTGGTGACGACGTTGACGTATGGCGTGCAGTGTTTTGCGCGTTATCCGGACCACATGAAACAACACGATTTCTTCAAATCTGCGATGCCGGAGGGTTACGTCCAGGAGCGTACCATTTCCTTCAAGGATGATGGCACTTACAAAACTCGCGCAGAGGTTAAGTTTGAAGGTGACACGCTGGTCAATCGTATCGAATTGAAGGGTATCGACTTTAAAGAGGATGGTAACATTCTGGGCCATAAACTGGAGTATAACTTCAACAGCCATAATGTTTACATTACGGCAGACAAGCAAAAGAACGGCATCAAGGCCAATTTCAAGATTCGCCACAATGTTGAGGACGGTAGCGTCCAACTGGCCGACCATTACCAGCAGAACACCCCAATTGGTGACGGTCCGGTTTTGCTGCCGGATAATCACTATCTGAGCACCCAAAGCGTGCTGAGCAAAGATCCGAACGAAAAACGTGATCACATGGTCCTGCTGGAATTTGTGACCGCTGCGGGCATCACCCACGGTATGGACGAGCTGTATAAGCGTCCGTAA", "sfGFP", "orf"
    ))
    parts.append(QuickBasicPart(
        "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCTTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA", "mCherry", "orf"
    ))
    parts.append(QuickBasicPart(
        "ATGTCCGAGTTGATCAAAGAGAACATGCATATGAAATTATATATGGAAGGCACTGTAGATAATCATCATTTTAAATGTACGTCGGAAGGCGAAGGTAAACCATATGAAGGTACGCAGACGATGCGCATCAAGGTGGTGGAGGGCGGTCCGCTGCCATTCGCTTTCGATATTTTAGCCACGAGCTTCCTCTACGGTTCTAAAACTTTCATCAATCACACGCAGGGTATTCCGGACTTCTTTAAACAGTCGTTCCCGGAGGGTTTCACCTGGGAACGCGTTACCACGTATGAAGATGGTGGTGTGCTTACGGCAACGCAGGACACGAGCCTTCAGGATGGGTGTTTGATTTACAACGTGAAAATTCGTGGTGTGAACTTCACGTCTAACGGCCCGGTGATGCAGAAAAAAACACTGGGTTGGGAAGCCTTTACCGAAACCCTGTATCCGGCGGACGGTGGCCTGGAAGGCCGTAATGATATGGCCTTGAAATTAGTCGGCGGTTCACACCTGATCGCGAACGCGAAAACAACCTATCGTAGTAAAAAACCAGCCAAAAACCTGAAAATGCCGGGCGTCTACTACGTAGACTACCGTCTGGAGCGCATTAAAGAGGCGAATAATGAAACCTATGTCGAGCAGCACGAAGTTGCGGTTGCACGCTATTGCGATCTGCCCAGCAAACTGGGCCACAAGCTTAATGGTAGCTAA", "mTagBFP2", "orf"
    ))
    for part in parts:
        make_gb_from_quick_part(part)


def add_feature_w_seq(seqrec, feature_seq, feature_type="misc_feature", **qualifiers):
    """
    Add a SeqFeature to a SeqRecord (seqrec) using that feature's sequence (feature_seq). Forward strand only.
    Each kwarg value in qualifiers must be given as a list e.g. key=["value"]

    """
    qualifier_dict = {qualifier[0]: qualifier[1] for qualifier in qualifiers.items()}
    feature_locations = nt_search(str(seqrec.seq), feature_seq)
    for location in feature_locations[1:]:
        seqrec.features.append(
            SeqFeature(
                type=feature_type,
                location=FeatureLocation(location, location + len(feature_seq), strand=+1),
                qualifiers=qualifier_dict
            )
        )


if __name__ == "__main__":
    main()