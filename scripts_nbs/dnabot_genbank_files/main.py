import basicsynbio as bsb
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqUtils import nt_search
from dataclasses import dataclass


def main():
    generate_vectors()


def generate_vectors():
    paths = [
        "basic_seva_18_ampr-puc-1.gb",
        "basic_seva_37_cmr-p15a-1.gb"
    ]
    vectors = (bsb.import_part(handle, "genbank") for handle in paths)
    linkers_features = (bsb.biolegio_dict["LMP"], bsb.biolegio_dict["LMS"])
    for vector in vectors:
        vector.id = f"{vector.id}1"
        vector.description="Vector for BASIC DNA assembly containing \
            mScarlet counter selection marker"
        vector.name = vector.id[:len(vector.id)-2]
        vector.annotations = bsb.DEFAULT_ANNOTATIONS
        for linker in linkers_features:
            add_feature_w_seq(
                vector,
                str(linker.features[0].extract(linker.seq)),
                function=["BASIC DNA assembly linker"],
                standard_name=[str(linker.id)],
                note=[str(linker.id)]
            )
        SeqIO.write(vector, f"{vector.id}.gb", "genbank")


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