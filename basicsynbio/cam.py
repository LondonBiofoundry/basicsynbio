"""Module contains a collection of objects for computer assisted manufacturing
within the BASIC DNA assembly framework."""

from Bio.SeqUtils import molecular_weight
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .main import (
    BasicAssembly,
    BasicLinker,
    BasicPart,
    BasicUTRRBSLinker,
    ClipReaction,
    LinkerException,
)
import csv
from dataclasses import dataclass
from datetime import datetime
from collections import OrderedDict, Counter
import hashlib
import json
import os
from pathlib import Path
from typing import Union, Dict, Generator
import re
import zipfile
import string
import math

WELLS_384 = [
    letter + str(number)
    for number in range(1, 25)
    for letter in string.ascii_uppercase[0:16]
]
WELLS_384_8_CHANNEL = [
    letter + str(number)
    for number in range(1, 25)
    for letter in string.ascii_uppercase[0:16:2]
]
WELLS_96 = [
    letter + str(number)
    for number in range(1, 13)
    for letter in string.ascii_uppercase[:8]
]


def new_part_resuspension(
    part: BasicPart, mass: float, double_stranded: bool = True
) -> float:
    """Returns the volume of resuspension buffer (µL) required for a 75 nM
    solution of part, equivalent to 75 fmol/µL.

    Args:
        part: BasicPart object.
        mass: mass of synthesised part (ng).
        double_stranded (optional): True (default) indicates part is dsDNA.
    """
    return (
        (mass * 10 ** -9)
        / molecular_weight(part.seq, double_stranded=double_stranded)
        * 1
        / (75e-9)
        * 10 ** 6
    )


class BasicBuild:
    """Class provides methods and attributes for building BasicAssembly
    objects.

    Attributes:
        basic_assemblies: A collection of all the
            assemblies within the BasicBuild object.
        clips_data: A dictionary of each clip reaction as a key
            with associated assemblies as values.
        unique_clips: A list of each unique ClipReaction
            within the BasicBuild object.
        unique_parts: Tuple of the unique BasicPart object required
            to implement the build.
        unique_linkers: Tuple of the unique BasicLinker objects
            required to implement the build.
    """

    def __init__(self, *basic_assemblies: BasicAssembly):
        """Initiate BasicBuild.

        Args:
            *basic_assemblies: BasicAssembly objects.
        """
        self.basic_assemblies = basic_assemblies
        self.clips_data = self._return_clips_data()
        self.unique_clips = [clip for clip in self.clips_data.keys()]
        self.unique_parts_data = self._unique_parts_linkers_data(
            "part", *(clip_reaction._part for clip_reaction in self.clips_data)
        )
        self.unique_linkers_data = self._unique_parts_linkers_data(
            "linker", *(clip_reaction._prefix for clip_reaction in self.clips_data)
        )
        for clip_reaction in self.unique_clips:
            part_hash = seqrecord_hexdigest(clip_reaction._part)
            prefix_hash = seqrecord_hexdigest(clip_reaction._prefix)
            self.unique_parts_data[part_hash]["clip_reactions"].append(clip_reaction)
            self.unique_linkers_data[prefix_hash]["clip_reactions"].append(
                clip_reaction
            )
        self.unique_parts = tuple(
            part_dict["part"] for part_dict in self.unique_parts_data.values()
        )
        self.unique_linkers = tuple(
            linker_dict["linker"] for linker_dict in self.unique_linkers_data.values()
        )

    def update_parts(self, *parts: BasicPart) -> None:
        """Updates BasicBuild instance with parts.

        Replaces all existing BasicParts used in assemblies with the matching equivalent in parts. Parts are considered equivalent if seqrecord_hexdigest() returns identical values for both.

        Args:
            parts: parts to replace the BasicParts used in
                assemblies

        Raises:
            ValueError: if length parts is not equal to lenght this build's
                unique parts
        """
        if len(parts) != len(self.unique_parts_data):
            raise ValueError(
                f"length of *parts is {len(parts)} whereas self.unqiue_parts has {len(self.unique_parts_data)} elements. The two must match."
            )
        parts_dict = self._unique_parts_linkers_data("part", *parts)
        basic_assemblies = []
        for assembly in self.basic_assemblies:
            parts_linkers = [
                part_linker
                if isinstance(part_linker, BasicLinker)
                else parts_dict[seqrecord_hexdigest(part_linker)]["part"]
                for part_linker in assembly.parts_linkers
            ]
            basic_assemblies.append(BasicAssembly(assembly.id, *parts_linkers))
        self.__init__(*basic_assemblies)

    def _return_clips_data(self) -> dict:
        """Analyses the build and returns a dictionary describing ClipReactions and their associated assemblies.

        Returns:
            clips_dict: each ClipReaction as a key with values describing
            associated assemblies.

        """
        clips_dict = OrderedDict(
            **{
                clip_reaction: []
                for assembly in self.basic_assemblies
                for clip_reaction in assembly.clip_reactions
            }
        )
        for assembly in self.basic_assemblies:
            for clip_reaction in assembly.clip_reactions:
                clips_dict[clip_reaction].append(assembly)
        return clips_dict

    def _unique_parts_linkers_data(
        self, object_key: str, *parts_linkers: Union[BasicPart, BasicLinker]
    ) -> Dict[str, dict]:
        """Returns a dictionary of unique objects in parts_linkers. Includes
        an empty list for each item to populate with clip_reactions used by
        each unique part/linker.

        Args:
            object_key: can be "part" or "linker".
            parts_linkers: collection of parts and linkers present within build
                used to search through for uniqueness.

        Returns:
            dict: hexdigest of part/linker along with a sub-dictionary containing part/linker and empty list for populating associated clip_reactions.
        """
        return {
            seqrecord_hexdigest(part_linker): {
                object_key: part_linker,
                "clip_reactions": [],
            }
            for part_linker in parts_linkers
        }

    def _duplicate_assembly_ids(self, assemblies):
        assemblies_ids = [assembly.id for assembly in assemblies]
        if len(set(assemblies_ids)) < len(assemblies):
            top_assembly_id = Counter(assemblies_ids).most_common(1)[0]
            raise BuildException(
                f"ID '{top_assembly_id[0]}' has been assigned to {top_assembly_id[1]} BasicAssembly instance/s. All assemblies of a build should have a unique 'id' attribute."
            )

    def export_csvs(self, path: str = None):
        """Writes information about each clip_data and assembly to
        two dependent CSV files in the same folder the command
        is executed.

        Args:
            path (optional): path to zipped folder of csv files. If none defaults to
                working directory with a time stamped name, output csvs is created.
        """
        if path == None:
            now = datetime.now()
            zip_path = Path.cwd() / f"build_{now.strftime('%d-%m-%Y_%H.%M.%S')}.zip"
        else:
            zip_path = path
        with open(Path.cwd() / "clips.csv", "w", newline="") as f:
            fieldnames = [
                "Clip Index",
                "Prefix ID",
                "Part ID",
                "Part Name",
                "Part suggested stock concentration (ng/µL)",
                "Part stock per 30 µL clip (µL)",
                "Suffix ID",
                "Total assemblies",
                "Assembly indexes",
            ]
            thewriter = csv.DictWriter(f, fieldnames=fieldnames)
            thewriter.writeheader()
            for index, clip_data in enumerate(self.clips_data.items()):
                thewriter.writerow(
                    {
                        "Clip Index": index + 1,
                        "Prefix ID": clip_data[0]._prefix.prefix_id,
                        "Part ID": clip_data[0]._part.id,
                        "Part Name": clip_data[0]._part.name,
                        "Part suggested stock concentration (ng/µL)": clip_data[
                            0
                        ]._part.concentration(),
                        "Part stock per 30 µL clip (µL)": 1,
                        "Suffix ID": clip_data[0]._suffix.suffix_id,
                        "Total assemblies": len(clip_data[1]),
                        "Assembly indexes": [
                            self.basic_assemblies.index(assembly) + 1
                            for assembly in clip_data[1]
                        ],
                    }
                )
        with open(Path.cwd() / "assemblies.csv", "w", newline="") as f:
            fieldnames = ["Assembly Index", "Assembly ID", "Clip indexes"]
            thewriter = csv.DictWriter(f, fieldnames=fieldnames)
            thewriter.writeheader()
            for index, assembly in enumerate(self.basic_assemblies):
                thewriter.writerow(
                    {
                        "Assembly Index": index + 1,
                        "Assembly ID": assembly.id,
                        "Clip indexes": [
                            self.unique_clips.index(clip_reaction) + 1
                            for clip_reaction in assembly.clip_reactions
                        ],
                    }
                )
        with zipfile.ZipFile(zip_path, "w") as my_zip:
            my_zip.write("assemblies.csv")
            my_zip.write("clips.csv")
        os.remove(Path.cwd() / "assemblies.csv")
        os.remove(Path.cwd() / "clips.csv")

    def export_echo_assembly(
        self,
        path: str = None,
        bufferWell: str = "A1",
        waterWell: str = "B1",
        alternate_well: bool = False,
    ) -> None:
        """Writes automation scripts for a echo liquid handler to build assemblies from clips.

        Args:
            path (optional): path to zipped folder of csv files. If none defaults to
                working directory with a time stamped name, output csvs is created.
            bufferWell (optional): location in 6 well plate of assembly buffer.
            waterWell (optional): location in 6 well plate of dH20.
            alternate_well (optional): specifies whether alternating wells are to be used in the input 384 well plate.

        Returns:
            str: Path of zip file containing echo automation scripts

        Raises:
            ValueError: If waterWell or BufferWell is not in ["A1", "B1", "A2", "B2", "A3", "B3"]; if self contains
                96 or more assemblies or if the build requires equal or more than 384 used clip wells for alternate_well(True)
                or 192 for alternate_well(False).
        """
        # Volumes given in nanoliters
        CLIP_VOLUME_FOLLOWING_MAGBEAD = 40000
        CLIP_VOLUME_PER_ASSEMBLY = 500
        BUFFER_VOLUME_PER_ASSEMBLY = 500
        TOTAL_ASSEMBLY_VOLUME = 5000
        DEADSPACE_96_WELL_PLATE = 20000
        CLIPS_PER_WELL = (
            CLIP_VOLUME_FOLLOWING_MAGBEAD - DEADSPACE_96_WELL_PLATE
        ) / CLIP_VOLUME_PER_ASSEMBLY
        if waterWell not in ["A1", "B1", "A2", "B2", "A3", "B3"]:
            raise ValueError(
                "Water Well location needs to be within the 6 well plate, between A1 - B3"
            )
        if bufferWell not in ["A1", "B1", "A2", "B2", "A3", "B3"]:
            raise ValueError(
                "Assembly Buffer Well location needs to be within the 6 well plate, between A1 - B3"
            )
        if len(self.basic_assemblies) >= 96:
            raise ValueError(
                """To many assemblies in the build to be handled by a single 96 
                well plate, use less than 96 assemblies within the build"""
            )

        well_384_index = (
            0  # index of WELLS_384(_8_CHANNEL) which increments as wells are assigned
        )
        # to clips facilitates multiple wells being assigned to one clip, as increments on new clip
        # and when chip is used more than each multiple of CLIPS_PER_WELL
        index_origin_plate_mapping = (
            {}
        )  # A dictionary with keys as index in clipdata object and
        # values of a dictionary containing key of wells and values of number of times used e.g.
        # > {0: {'A1': 0, 'B1': 0}, 1: {'C1': 0},... This says clip at index 0 is in wells A1 & B1.
        # > {0: {'A1': 40, 'B1': 23}, 1: {'C1': 1}... After the CSV writing loop the dict looks
        # similar to this as for clip index 0 A1 is used until is has been used CLIPS_PER_WELL times
        # then the next index is used in this case B1.

        # Mapping clips to associated 384 plate wells considing wells needed from volume
        # and alternate_well attribute
        for index, clip_data in enumerate(self.clips_data.items()):
            wells_for_current_clip = {}
            for i in range(1 + math.floor(len(clip_data[1]) / (CLIPS_PER_WELL))):
                if not alternate_well:
                    wells_for_current_clip.update({WELLS_384[well_384_index]: 0})
                else:
                    wells_for_current_clip.update(
                        {WELLS_384_8_CHANNEL[well_384_index]: 0}
                    )
                well_384_index += 1
                if well_384_index >= (384 - alternate_well * 192):
                    raise ValueError(
                        """The Number of clip wells used cannot exceed 384 for alternate_well(True) or 192 for
                        alternate_well(False), please choose a build with less unique clips or change the
                        alternate_well setting"""
                    )
            index_origin_plate_mapping.update({index: wells_for_current_clip})

        if path == None:
            now = datetime.now()
            zip_path = (
                Path.cwd()
                / f"Echo_Instructions_{now.strftime('%d-%m-%Y_%H.%M.%S')}.zip"
            )
        else:
            zip_path = path
        # Create list of transfers by connecting clips used in assembly to clips wells
        # ensuring each well is used maximum CLIPS_PER_WELL times
        with open(Path.cwd() / "echo_clips_1.csv", "w", newline="") as f1, open(
            Path.cwd() / "echo_water_buffer_2.csv", "w", newline=""
        ) as f2:
            fieldnames = ["Destination Well", "Source Well", "Transfer Volume"]
            thewriter_clips = csv.DictWriter(f1, fieldnames=fieldnames)
            thewriter_clips.writeheader()
            thewriter_water_buffer = csv.DictWriter(f2, fieldnames=fieldnames)
            thewriter_water_buffer.writeheader()
            for index, assembly in enumerate(self.basic_assemblies):
                for clip in [
                    self.unique_clips.index(clip_reaction)
                    for clip_reaction in assembly.clip_reactions
                ]:
                    for well in index_origin_plate_mapping[clip].keys():
                        if index_origin_plate_mapping[clip][well] < CLIPS_PER_WELL:
                            index_origin_plate_mapping[clip][well] += 1
                            break
                    thewriter_clips.writerow(
                        {
                            "Destination Well": WELLS_96[index],
                            "Source Well": well,
                            "Transfer Volume": CLIP_VOLUME_PER_ASSEMBLY,
                        }
                    )
                thewriter_water_buffer.writerow(
                    {
                        "Destination Well": WELLS_96[index],
                        "Source Well": bufferWell,
                        "Transfer Volume": BUFFER_VOLUME_PER_ASSEMBLY,
                    }
                )
                thewriter_water_buffer.writerow(
                    {
                        "Destination Well": WELLS_96[index],
                        "Source Well": waterWell,
                        "Transfer Volume": TOTAL_ASSEMBLY_VOLUME
                        - BUFFER_VOLUME_PER_ASSEMBLY
                        - CLIP_VOLUME_PER_ASSEMBLY
                        * len(
                            [
                                self.unique_clips.index(clip_reaction)
                                for clip_reaction in assembly.clip_reactions
                            ]
                        ),
                    }
                )
        index_origin_plate_mapping = {}
        with zipfile.ZipFile(zip_path, "w") as my_zip:
            my_zip.write("echo_clips_1.csv")
            my_zip.write("echo_water_buffer_2.csv")
        os.remove(Path.cwd() / "echo_clips_1.csv")
        os.remove(Path.cwd() / "echo_water_buffer_2.csv")
        return zip_path

    @property
    def basic_assemblies(self):
        return self._basic_assemblies

    @basic_assemblies.setter
    def basic_assemblies(self, values):
        if not all(issubclass(type(value), BasicAssembly) for value in values):
            raise TypeError("Not all *basic_assemblies are BasicAssembly instances.")
        self._duplicate_assembly_ids(values)
        self._basic_assemblies = values


class BuildEncoder(json.JSONEncoder):
    """A Class to encode BasicBuild objects extending `json.JSONEncoder` class"""

    def default(self, obj):
        if isinstance(obj, BasicBuild):
            return {
                "unique_parts": self.unique_parts_json(obj),
                "unique_linkers": self.unique_linkers_json(obj),
                "clips_data": self.clips_data_json(obj),
                "assembly_data": self.assembly_data_json(obj),
                "__BasicBuild__": True,
            }
        return super().default(obj)

    @staticmethod
    def unique_parts_json(obj):
        """A function to create machine readable json objects, completely
        describing unique parts.

        Args:
            obj: BasicBuild object to be decoded

        Returns:
            dictionary: Each item describes a unique BasicPart object required for the build.
        """
        return {
            key: {
                "sequence": str(value["part"].seq),
                "id": value["part"].id,
                "name": value["part"].name,
                "description": value["part"].description,
                "suggested stock concentration (ng/µL)": value["part"].concentration(),
                "stock per 30 µL clip (µL)": 1,
                "total clip reactions": len(value["clip_reactions"]),
                "clip reactions": [
                    clip_reaction._hexdigest()
                    for clip_reaction in value["clip_reactions"]
                ],
            }
            for key, value in obj.unique_parts_data.items()
        }

    @staticmethod
    def unique_linkers_json(obj):
        """A function to create machine readable json objects, completely
        describing unique linkers.

        Args:
            obj: BasicBuild object to be decoded

        Returns:
            dictionary: Each item describes a given unique BasicLinker object required for the build.
        """
        return {
            key: {
                "id": value["linker"].id,
                "linker_class": str(type(value["linker"])),
                "sequence": str(value["linker"].seq),
                "prefix_id": value["linker"].prefix_id,
                "suffix_id": value["linker"].suffix_id,
                "total clip reactions": len(value["clip_reactions"]),
                "clip reactions": [
                    clip_reaction._hexdigest()
                    for clip_reaction in value["clip_reactions"]
                ],
            }
            for key, value in obj.unique_linkers_data.items()
        }

    @staticmethod
    def clips_data_json(obj):
        """A function to create machine readable json objects, completely
        describing ClipReaction objects within BasicBuild.

        Args:
            obj: BasicBuild object to be decoded

        Returns:
            dictionary: Each item describes a given unique ClipReaction object required for the build.
        """
        return {
            key._hexdigest(): {
                "prefix": {
                    "key": seqrecord_hexdigest(key._prefix),
                    "prefix_id": key._prefix.prefix_id,
                },
                "part": {
                    "key": seqrecord_hexdigest(key._part),
                    "id": key._part.id,
                    "name": key._part.name,
                },
                "suffix": {
                    "key": seqrecord_hexdigest(key._suffix),
                    "suffix_id": key._suffix.suffix_id,
                },
                "total assemblies": len(value),
                "assembly indexes": [
                    obj.basic_assemblies.index(assembly) for assembly in value
                ],
            }
            for key, value in obj.clips_data.items()
        }

    @staticmethod
    def assembly_data_json(obj):
        """A function to create machine readable json objects, completely
        describing BasicAssembly objects within BasicBuild.

        Args:
            obj: BasicBuild object to be decoded

        Returns:
            list: returns a list of dictionaries populated by
            dictionaries with items describing BasicAssemblies generated by the build.
        """
        return [
            {
                "id": assembly.id,
                "clip reactions": [
                    clip_reaction._hexdigest()
                    for clip_reaction in assembly.clip_reactions
                ],
            }
            for assembly in obj.basic_assemblies
        ]


class BuildDecoder(json.JSONDecoder):
    """A Class to decode json dictionary to basicsynbio objects,
    extending `json.JSONDecoder` class

    """

    def __init__(self):
        json.JSONDecoder.__init__(self, object_hook=self.decode_build)

    def decode_build(self, dictionary: dict) -> BasicBuild:
        """A method to return BasicBuild from encoded json object

        Args:
            dictionary: json object encoded by BuildEncoder

        Returns:
            BasicBuild: BasicBuild object built from encoded json
        """
        if "__BasicBuild__" in dictionary:
            self.unique_parts_data = self.return_unqiue_parts(dictionary)
            self.unique_linkers_data = self.return_unique_linkers(dictionary)
            basic_assemblies = self.return_basic_assemblies(dictionary)
            return BasicBuild(*basic_assemblies)
        return dictionary

    def return_basic_assemblies(
        self, dictionary: dict
    ) -> Generator[BasicAssembly, None, None]:
        """A method to yield BasicAssembly objects from encoded json object

        Args:
            dictionary: json object encoded by BuildEncoder

        Yields:
            BasicAssembly: BasicAssembly objects within encoded BasicBuild
        """
        for assembly in dictionary["assembly_data"]:
            parts_linkers = []
            for clip_reaction in assembly["clip reactions"]:
                parts_linkers += [
                    self.unique_linkers_data[
                        dictionary["clips_data"][clip_reaction]["prefix"]["key"]
                    ],
                    self.unique_parts_data[
                        dictionary["clips_data"][clip_reaction]["part"]["key"]
                    ],
                ]
            yield BasicAssembly(assembly["id"], *parts_linkers)

    @staticmethod
    def return_unqiue_parts(dictionary):
        """A method to return unique BasicPart objects from encoded json object

        Args:
            dictionary: json object encoded by BuildEncoder

        Returns:
            dictionary: containing unique BasicPart objects
        """
        return {
            key: BasicPart(
                seq=Seq(value["sequence"]),
                id=value["id"],
                name=value["name"],
                description=value["description"],
            )
            for key, value in dictionary["unique_parts"].items()
        }

    @staticmethod
    def return_unique_linkers(dictionary):
        """A method to return unqiue BasicLinker objects from encoded json object

        Args:
            dictionary: json object encoded by BuildEncoder

        Returns:
            dictionary: containing unique BasicLinker objects
        """
        unique_linkers = {}
        for key, value in dictionary["unique_linkers"].items():
            if re.match(".*BasicLinker", value["linker_class"]):
                unique_linker = BasicLinker(seq=Seq(value["sequence"]), id=value["id"])
            elif re.match(".*BasicUTRRBSLinker", value["linker_class"]):
                unique_linker = BasicUTRRBSLinker(
                    seq=Seq(value["sequence"]), id=value["id"]
                )
            else:
                raise LinkerException(
                    f"unique linker '{key}' does not have a recognised 'linker_class' attribute."
                )
            unique_linker.prefix_id = value["prefix_id"]
            unique_linker.suffix_id = value["suffix_id"]
            unique_linkers[key] = unique_linker
        return unique_linkers


class BuildException(Exception):
    pass


def seqrecord_hexdigest(seqrecord_obj: Union[SeqRecord, BasicPart, BasicLinker]) -> str:
    """Returns an MD5 hash of a Bio.SeqRecord.SeqRecord-like object, using
    relevant attributes.

    Note:
        Care should be taken when using the return value of this function as an identifier/hash for seqrecord_obj given these objects are mutable.

    Args:
        seqrecord_obj: Bio.SeqRecord.SeqRecord-like object, containing relevant
            attributes

    Returns:
        MD5 hash
    """
    seqrec_hash = hashlib.md5(str(seqrecord_obj.seq).encode("UTF-8"))
    bytes_objs = [
        getattr(seqrecord_obj, attribute).encode("UTF-8")
        for attribute in ["name", "description"]
    ]
    for element in bytes_objs:
        seqrec_hash.update(element)
    return seqrec_hash.hexdigest()
