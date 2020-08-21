"""Module contains a collection of objects for computer assisted manufacturing within the BASIC DNA assembly framework."""

from Bio.SeqUtils import molecular_weight
from Bio.Seq import Seq
from .main import (
    BasicAssembly,
    BasicLinker,
    BasicPart,
    BasicUTRRBSLinker,
    ClipReaction,
    LinkerException    
)
from dataclasses import dataclass
from collections import OrderedDict, Counter
import json
import hashlib
import re


def new_part_resuspension(part, mass: float, double_stranded=True):
    """Returns the volume of resuspension buffer (µL) required for a 75 nM solution of part, equivalent to 75 fmol/µL.

    Args:
        part -- BasicPart object.
        mass -- mass of synthesised part (ng).
        double_stranded -- True (default) indicates part is dsDNA.

    """
    return (mass*10**-9)/molecular_weight(part.seq, double_stranded=double_stranded)*1/(75e-9)*10**6


class BasicBuild():
    """Class provides methods and attributes for building BasicAssembly objects."""

    def __init__(self, *basic_assemblies):
        """Initiate BasicBuild.

        Args:
            *basic_assemblies -- BasicAssembly objects.

        """
        self.basic_assemblies = basic_assemblies
        self.clips_data = self._return_clips_data()
        self.unique_parts = self._unique_parts_linkers(
            "part",
            *(clip_reaction._part for clip_reaction in self.clips_data)
        )
        self.unique_linkers = self._unique_parts_linkers(
            "linker",
            *(clip_reaction._prefix for clip_reaction in self.clips_data)
        )
        for clip_reaction in self.clips_data.keys():
            part_hash = _seqrecord_hexdigest(clip_reaction._part)
            prefix_hash = _seqrecord_hexdigest(clip_reaction._prefix)
            self.unique_parts[part_hash]["clip_reactions"].append(clip_reaction)
            self.unique_linkers[prefix_hash]["clip_reactions"].append(clip_reaction)

    def _return_clips_data(self):
        """Returns a dictionary of ClipReactions with values describing basic_assemblies it uses."""
        clips_dict = OrderedDict(
            **{clip_reaction: [] for assembly in self.basic_assemblies for clip_reaction in assembly.clip_reactions})
        for assembly in self.basic_assemblies:
            for clip_reaction in assembly.clip_reactions:
                clips_dict[clip_reaction].append(assembly)
        return clips_dict

    def _unique_parts_linkers(self, object_key: str, *parts_linkers):
        """Returns a dictionary of unique objects based on hash of 'id' and 'seq' attribute. Includes an empty list associated to populate with clip_reactions used by each unique part/linker.
        
        Args:
            object_key -- "part" or "linker".
        """
        return {
            _seqrecord_hexdigest(part_linker): {
                object_key: part_linker,
                "clip_reactions": []
            } for part_linker in parts_linkers
        }

    def _duplicate_assembly_ids(self, assemblies):
        """If multiple elements of self.basic_assemblies have same "id" attribute, raises a BuildException"""
        assemblies_ids = [assembly.id for assembly in assemblies]
        if len(set(assemblies_ids)) < len(assemblies):
            top_assembly_id = Counter(assemblies_ids).most_common(1)[0]
            raise BuildException(
                f"ID '{top_assembly_id[0]}' has been assigned to {top_assembly_id[1]} BasicAssembly instance/s. All assemblies of a build should have a unique 'id' attribute.")

    @property
    def basic_assemblies(self):
        return self._basic_assemblies

    @basic_assemblies.setter
    def basic_assemblies(self, values):
        if not all(
                issubclass(type(value), BasicAssembly) for value in values):
            raise TypeError(
                "Not all *basic_assemblies are BasicAssembly instances."
            )
        self._duplicate_assembly_ids(values)
        self._basic_assemblies = values


class BuildEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, BasicBuild):
            return {
                "unique_parts": self.unique_parts_json(obj),
                "unique_linkers": self.unique_linkers_json(obj),
                "clips_data": self.clips_data_json(obj),
                "assembly_data": self.assembly_data_json(obj),
                "__BasicBuild__": True
            }
        return super().default(obj)
    
    @staticmethod
    def unique_parts_json(obj):
        return {
            key: {
                "sequence": str(value["part"].seq),
                "id": value["part"].id,
                "name": value["part"].name,
                "description": value["part"].description,
                "clip_reactions": [clip_reaction._hexdigest() for clip_reaction in value["clip_reactions"]]
            }
        for key, value in obj.unique_parts.items()}
    
    @staticmethod
    def unique_linkers_json(obj):
        return {
            key: {
                "id": value["linker"].id,
                "linker_class": str(type(value["linker"])),
                "sequence": str(value["linker"].seq),
                "prefix_id": value["linker"].prefix_id,
                "suffix_id": value["linker"].suffix_id,
                "clip_reactions": [clip_reaction._hexdigest() for clip_reaction in value["clip_reactions"]]
            }
        for key, value in obj.unique_linkers.items()}
    
    @staticmethod
    def clips_data_json(obj):
        return {
            key._hexdigest(): {
                "prefix": {
                    "key": _seqrecord_hexdigest(key._prefix),
                    "id": key._prefix.prefix_id
                },
                "part": {
                    "key": _seqrecord_hexdigest(key._part),
                    "id": key._part.id,
                    "name": key._part.name
                },
                "suffix": {
                    "key": _seqrecord_hexdigest(key._suffix),
                    "id": key._suffix.suffix_id
                },
                "assembly_data_indexes": [obj.basic_assemblies.index(assembly) for assembly in value]
            }
        for key, value in obj.clips_data.items()}
    
    @staticmethod
    def assembly_data_json(obj):
        return [
            {
                "id": assembly.id,
                "clip_reactions": [clip_reaction._hexdigest() for clip_reaction in assembly.clip_reactions]
            }
        for assembly in obj.basic_assemblies]


class BuildException(Exception):
    pass


class BuildDecoder(json.JSONDecoder):
    def __init__(self):
        json.JSONDecoder.__init__(self, object_hook=self.decode_build)

    def decode_build(self, dictionary):
        if "__BasicBuild__" in dictionary:
            self.unique_parts = self.return_unqiue_parts(dictionary)
            self.unique_linkers = self.return_unique_linkers(dictionary)
            basic_assemblies = self.return_basic_assemblies(dictionary)
            return BasicBuild(*basic_assemblies)
        return dictionary

    def return_basic_assemblies(self, dictionary):
        for assembly in dictionary["assembly_data"]:
            parts_linkers = []
            for clip_reaction in assembly["clip_reactions"]:
                parts_linkers += [
                    self.unique_linkers[dictionary["clips_data"][clip_reaction]["prefix"]["key"]],
                    self.unique_parts[dictionary["clips_data"][clip_reaction]["part"]["key"]],
                ]
            yield BasicAssembly(assembly["id"], *parts_linkers)

    @staticmethod
    def return_unqiue_parts(dictionary):
        return {
            key: BasicPart(
                seq=Seq(value["sequence"]),
                id=value["id"],
                name=value["name"],
                description=value["description"]
            )
        for key, value in dictionary["unique_parts"].items()}

    @staticmethod
    def return_unique_linkers(dictionary):
        unique_linkers = {}
        for key, value in dictionary["unique_linkers"].items():
            if re.match(".*BasicLinker", value["linker_class"]):
                unique_linker = BasicLinker(
                    seq=Seq(value["sequence"]),
                    id=value["id"]
                )
            elif re.match(".*BasicUTRRBSLinker", value["linker_class"]):
                unique_linker = BasicUTRRBSLinker(
                    seq=Seq(value["sequence"]),
                    id=value["id"]
                )
            else:
                raise LinkerException(f"unique linker '{key}' does not have a recognised 'linker_class' attribute.")
            unique_linker.prefix_id = value["prefix_id"]
            unique_linker.suffix_id = value["suffix_id"]
            unique_linkers[key] = unique_linker
        return unique_linkers


def _seqrecord_hexdigest(seqrecord_obj):
        """Returns an MD5 hash of a Bio.SeqRecord.SeqRecord-like object, using relevant attributes."""
        seqrec_hash = hashlib.md5(str(seqrecord_obj.seq).encode("UTF-8"))
        bytes_objs = [getattr(seqrecord_obj, attribute).encode("UTF-8") for attribute in [
            "id",
            "name",
            "description"
        ]]
        for element in bytes_objs:
            seqrec_hash.update(element)
        return seqrec_hash.hexdigest()
