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


@dataclass
class ClipInfo:
    clip_reaction: ClipReaction
    basic_assemblies: tuple


class BasicBuild():
    """Class provides methods and attributes for building BasicAssembly objects."""

    def __init__(self, *basic_assemblies):
        """Initiate BasicBuild.

        Args:
            *basic_assembiles -- BasicAssembly objects.

        """
        self.basic_assemblies = basic_assemblies
        self.clips_info = self.return_clips_info()
        self.unique_parts = tuple(
            element.clip_reaction._part for element in self.clips_info)
        self.unique_linkers = tuple(
            element.clip_reaction._prefix for element in self.clips_info)
        self.clip_indexes = {element.clip_reaction: ind for ind,
                           element in enumerate(self.clips_info)}

    def return_clips_info(self):
        """Returns {ClipReaction: [BasicAssemblies]} where BasicAssembly objects in list require ClipReaction."""
        clips_dict = OrderedDict(
            **{clip_reaction: [] for assembly in self.basic_assemblies for clip_reaction in assembly.clip_reactions})
        for assembly in self.basic_assemblies:
            for clip_reaction in assembly.clip_reactions:
                clips_dict[clip_reaction].append(assembly)
        return tuple(ClipInfo(key, tuple(value)) for key, value in clips_dict.items())

    def _duplicate_assembly_ids(self, assemblies):
        """If multiple elements of assemblies have same "id" attribute, raises a BuildException"""
        assemblies_ids = [assembly.id for assembly in assemblies]
        if len(set(assemblies_ids)) < len(assemblies):
            top_assembly_id = Counter(assemblies_ids).most_common(1)[0]
            raise BuildException(
                f"ID '{top_assembly_id[0]}' has been assigned to {top_assembly_id[1]} BasicAssembly instance/s.")
    
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
