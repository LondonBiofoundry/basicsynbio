"""Module contains a collection of objects for computer assisted manufacturing within the BASIC DNA assembly framework."""

from Bio.SeqUtils import molecular_weight
from .main import BasicAssembly, ClipReaction
from dataclasses import dataclass
from collections import OrderedDict, Counter
import json


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
            part_hash = _seqrecord_hash(clip_reaction._part)
            prefix_hash = _seqrecord_hash(clip_reaction._prefix)
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
            _seqrecord_hash(part_linker): {
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
                "assembly_data": self.assembly_data_json(obj)
            }
        return super().default(obj)
    
    def unique_parts_json(self, obj):
        return {
            str(key): {
                "id": value["part"].id,
                "sequence": str(value["part"].seq),
                "clip_reactions": [str(hash(clip_reaction)) for clip_reaction in value["clip_reactions"]]
            }
        for key, value in obj.unique_parts.items()}
    
    def unique_linkers_json(self, obj):
        return {
            str(key): {
                "id": value["linker"].id,
                "linker_class": str(type(value["linker"])),
                "sequence": str(value["linker"].seq),
                "prefix_id": value["linker"].prefix_id,
                "suffix_id": value["linker"].suffix_id,
                "clip_reactions": [str(hash(clip_reaction)) for clip_reaction in value["clip_reactions"]]
            }
        for key, value in obj.unique_linkers.items()}

    def clips_data_json(self, obj):
        return {
            str(hash(key)): {
                "prefix": str(_seqrecord_hash(key._prefix)),
                "part": str(_seqrecord_hash(key._part)),
                "suffix": str(_seqrecord_hash(key._suffix)),
                "assembly_data_indexes": [obj.basic_assemblies.index(assembly) for assembly in value]
            }
        for key, value in obj.clips_data.items()}
        
    def assembly_data_json(self, obj):
        return [
            {
                "id": assembly.id,
                "clip_reactions": [str(hash(clip_reaction)) for clip_reaction in assembly.clip_reactions]
            }
        for assembly in obj.basic_assemblies]


class BuildException(Exception):
    pass


def _seqrecord_hash(seqrecord_obj):
        """Returns a hash of a Bio.SeqRecord.SeqRecord-like object."""
        return hash((seqrecord_obj.id, seqrecord_obj.seq))
