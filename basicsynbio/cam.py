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
            *basic_assembiles -- BasicAssembly objects.

        """
        self.basic_assemblies = basic_assemblies
        self.clips_data = self.return_clips_data()
        self.unique_parts = self._unique_parts_linkers(
            "part",
            *(element["clip_reaction"]._part for element in self.clips_data)
        )
        self.unique_linkers = self._unique_parts_linkers(
            "linker",
            *(element["clip_reaction"]._prefix for element in self.clips_data)
        )
        for ind, element in enumerate(self.clips_data):
            clip_reaction = element["clip_reaction"]
            part_hash = self._part_linker_hash(clip_reaction._part)
            prefix_hash = self._part_linker_hash(clip_reaction._prefix)
            self.unique_parts[part_hash]["clips_indexes"].append(ind)
            self.unique_linkers[prefix_hash]["clips_indexes"].append(ind)

    def return_clips_data(self):
        """Returns {ClipReaction: [BasicAssemblies]} where BasicAssembly objects in list require ClipReaction."""
        clips_dict = OrderedDict(
            **{clip_reaction: [] for assembly in self.basic_assemblies for clip_reaction in assembly.clip_reactions})
        for assembly in self.basic_assemblies:
            for clip_reaction in assembly.clip_reactions:
                clips_dict[clip_reaction].append(assembly)
        return tuple({"clip_reaction": key, "basic_assemblies": value} for key, value in clips_dict.items())

    def _unique_parts_linkers(self, object_key: str, *parts_linkers):
        """Returns a dictionary of unique objects based on hash of 'id' and 'seq' attribute. Includes an empty list associated to populate with clips_indexes.
        
        Args:
            object_key -- "part" or "linker".
        """
        return {
            self._part_linker_hash(part_linker): {
                object_key: part_linker,
                "clips_indexes": []
            } for part_linker in parts_linkers
        }

    def _part_linker_hash(self, part_linker):
        """Returns a hash of the part_linker based on id and sequence."""
        return hash((part_linker.id, part_linker.seq))

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


class BuildException(Exception):
    pass
