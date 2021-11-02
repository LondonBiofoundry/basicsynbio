"""Module contains key objects for computer assisted manufacturing (CAM)
within the BASIC DNA assembly framework.
"""

from Bio.Seq import Seq
from basicsynbio.main import (
    BasicAssembly,
    BasicLinker,
    BasicPart,
    BasicUTRRBSLinker,
    LinkerException,
)
from collections import OrderedDict, Counter
import json
from typing import Dict, Generator
import re


class BasicBuild:
    """Class provides methods and attributes for building BasicAssembly
    objects.

    Attributes:
        basic_assemblies: Tuple of BasicAssembly objects used to initiate BasicBuild.
        unique_clips: Tuple of unique ClipReaction objects
            required to implement the build.
        clips_data: A dictionary with each ClipReaction as a key
            and BasicAssembly objects that require it as values.
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

    def update_parts(self, *parts: BasicPart) -> None:
        """Updates BasicBuild instance with parts.

        Replaces all existing BasicParts used in assemblies with the matching equivalent in parts.

        Args:
            parts: parts to replace the BasicParts used in
                assemblies

        Raises:
            ValueError: if length parts is not equal to length this build's
                unique parts
        """
        if len(parts) != len(self.unique_parts_data):
            raise ValueError(
                f"length of *parts is {len(parts)} whereas self.unqiue_parts has {len(self.unique_parts_data)} elements. The two must match."
            )
        parts_dict = self._unique_parts_data(*parts)
        basic_assemblies = []
        for assembly in self.basic_assemblies:
            parts_linkers = [
                part_linker
                if isinstance(part_linker, BasicLinker)
                else parts_dict[part_linker.seq]["part"]
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
                for clip_reaction in assembly._clip_reactions
            }
        )
        for assembly in self.basic_assemblies:
            for clip_reaction in assembly._clip_reactions:
                clips_dict[clip_reaction].append(assembly)
        return clips_dict

    def _unique_parts_data(self, *parts: BasicPart) -> Dict[Seq, dict]:
        """Returns a dictionary of unique objects in parts. Includes
        an empty list for each item to populate with clip_reactions used by
        each unique part.

        Args:
            parts: collection of parts present within build
                used to search through for uniqueness.

        Returns:
            dict: seq attribute of part along with a sub-dictionary containing part and empty list for populating associated clip_reactions.
        """
        return {
            part.seq: {
                "part": part,
                "clip reactions": [],
            }
            for part in parts
        }

    def _unique_linkers_data(self, *linkers: BasicLinker) -> Dict[Seq, dict]:
        """Returns a dictionary of unique objects in linkers. Includes
        an empty lists to populate with clip_reactions using prefix or suffix linkers.

        Args:
            linkers: collection of linkers present within build
                used to search through for uniqueness.

        Returns:
            dict: seq attribute of linker along with a sub-dictionary containing linker and empty list for populating associated clip_reactions.
        """
        return {
            linker.seq: {
                "linker": linker,
                "prefix clip reactions": [],
                "suffix clip reactions": [],
            }
            for linker in linkers
        }

    def _set_unique_parts_linkers_data(self) -> None:
        """sets unique_part_data and unique_linker_data attributes."""
        self.unique_parts_data = self._unique_parts_data(
            *(clip_reaction._part for clip_reaction in self.clips_data)
        )
        self.unique_linkers_data = self._unique_linkers_data(
            *(clip_reaction._prefix for clip_reaction in self.clips_data)
        )
        for clip_reaction in self.unique_clips:
            self.unique_parts_data[clip_reaction._part.seq]["clip reactions"].append(
                clip_reaction
            )
            self.unique_linkers_data[clip_reaction._prefix.seq][
                "prefix clip reactions"
            ].append(clip_reaction)
            self.unique_linkers_data[clip_reaction._suffix.seq][
                "suffix clip reactions"
            ].append(clip_reaction)

    def _duplicate_assembly_ids(self, assemblies):
        assemblies_ids = [assembly.id for assembly in assemblies]
        if len(set(assemblies_ids)) < len(assemblies):
            top_assembly_id = Counter(assemblies_ids).most_common(1)[0]
            raise BuildException(
                f"ID '{top_assembly_id[0]}' has been assigned to {top_assembly_id[1]} BasicAssembly instance/s. All assemblies of a build should have a unique 'id' attribute."
            )

    @property
    def basic_assemblies(self):
        return self._basic_assemblies

    @basic_assemblies.setter
    def basic_assemblies(self, values):
        if not all(issubclass(type(value), BasicAssembly) for value in values):
            raise TypeError("Not all *basic_assemblies are BasicAssembly instances.")
        self._duplicate_assembly_ids(values)
        self._basic_assemblies = values
        self.clips_data = self._return_clips_data()
        self.unique_clips = tuple(clip for clip in self.clips_data.keys())
        self._set_unique_parts_linkers_data()
        self.unique_parts = tuple(
            part_dict["part"] for part_dict in self.unique_parts_data.values()
        )
        self.unique_linkers = tuple(
            linker_dict["linker"] for linker_dict in self.unique_linkers_data.values()
        )


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
            "UP"
            + str(index): {
                "sequence": str(value["part"].seq),
                "id": value["part"].id,
                "name": value["part"].name,
                "description": value["part"].description,
                "Part mass for 30 μL clip reaction (ng)": value["part"].clip_mass(),
                "clip reactions": [
                    "CR" + str(list(obj.clips_data.keys()).index(clip_reaction))
                    for clip_reaction in value["clip reactions"]
                ],
            }
            for index, value in enumerate(obj.unique_parts_data.values())
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
            "UL"
            + str(index): {
                "id": value["linker"].id,
                "linker_class": str(type(value["linker"])),
                "name": value["linker"].name,
                "sequence": str(value["linker"].seq),
                "prefix_id": value["linker"].prefix_id,
                "suffix_id": value["linker"].suffix_id,
                "prefix clip reactions": [
                    "CR" + str(list(obj.clips_data.keys()).index(clip_reaction))
                    for clip_reaction in value["prefix clip reactions"]
                ],
                "suffix clip reactions": [
                    "CR" + str(list(obj.clips_data.keys()).index(clip_reaction))
                    for clip_reaction in value["suffix clip reactions"]
                ],
            }
            for index, value in enumerate(obj.unique_linkers_data.values())
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
        linker_seqs = [linker_seq for linker_seq in obj.unique_linkers_data.keys()]
        part_seqs = [part_seq for part_seq in obj.unique_parts_data.keys()]
        return {
            "CR"
            + str(index): {
                "prefix": {
                    "key": "UL" + str(linker_seqs.index(value[0]._prefix.seq)),
                    "prefix_id": value[0]._prefix.prefix_id,
                },
                "part": {
                    "key": "UP" + str(part_seqs.index(value[0]._part.seq)),
                    "id": value[0]._part.id,
                    "name": value[0]._part.name,
                },
                "suffix": {
                    "key": "UL" + str(linker_seqs.index(value[0]._suffix.seq)),
                    "suffix_id": value[0]._suffix.suffix_id,
                },
                "total assemblies": len(value[1]),
                "assembly keys": [
                    "A" + str(obj.basic_assemblies.index(assembly))
                    for assembly in value[1]
                ],
            }
            for index, value in enumerate(obj.clips_data.items())
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
        return {
            "A"
            + str(index): {
                "id": assembly.id,
                "clip reactions": [
                    "CR" + str(list(obj.clips_data.keys()).index(clip_reaction))
                    for clip_reaction in assembly._clip_reactions
                ],
            }
            for index, assembly in enumerate(obj.basic_assemblies)
        }


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
        for assembly in dictionary["assembly_data"].values():
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
                print("part/linkers", parts_linkers)
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
                unique_linker = BasicLinker(
                    seq=Seq(value["sequence"]),
                    id=value["id"],
                    name=value["name"],
                )
            elif re.match(".*BasicUTRRBSLinker", value["linker_class"]):
                unique_linker = BasicUTRRBSLinker(
                    seq=Seq(value["sequence"]),
                    id=value["id"],
                    name=value["name"],
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
