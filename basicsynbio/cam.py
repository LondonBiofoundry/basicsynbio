"""Module contains a collection of objects for computer assisted manufacturing within the BASIC DNA assembly framework."""

from Bio.SeqUtils import molecular_weight
from .main import BasicAssembly


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
        all_parts = []
        all_linkers = []
        for assembly in self.basic_assemblies:
            all_parts += [clip_reaction.part for clip_reaction in assembly.clip_reactions]
            all_linkers += [clip_reaction.prefix for clip_reaction in assembly.clip_reactions]
        self.unique_parts = list(self._unique_objects_by_id(*all_parts))
        self.unique_linkers = list(self._unique_objects_by_id(*all_linkers))
    
    def _unique_objects_by_id(self, *instances):
        """From *instances, returns a collection of elements which all have a unique id attribute."""
        all_ids = [instance.id for instance in instances]
        unqiue_ids = {instance.id for instance in instances}
        for unique_id in unqiue_ids:
            yield instances[all_ids.index(unique_id)]

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
        self._basic_assemblies = values
