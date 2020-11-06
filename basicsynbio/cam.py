"""Module contains a collection of objects for computer assisted manufacturing
within the BASIC DNA assembly framework."""

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
import csv


def new_part_resuspension(part, mass: float, double_stranded=True):
    """Returns the volume of resuspension buffer (µL) required for a 75 nM
    solution of part, equivalent to 75 fmol/µL.

    Args:
        part -- BasicPart object.
        mass -- mass of synthesised part (ng).
        double_stranded -- True (default) indicates part is dsDNA.
    """
    return (mass*10**-9)/molecular_weight(part.seq, double_stranded=double_stranded)*1/(75e-9)*10**6


class BasicBuild():
    """Class provides methods and attributes for building BasicAssembly
    objects."""

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
    
    def update_parts(self, *parts):
        """Updates BasicBuild instance with *parts, replacing existing all
        BasicParts used in assemblies with the matching equivalent in.

        *parts.
        """
        if len(parts) != len(self.unique_parts):
            raise ValueError(f"length of *parts is {len(parts)} whereas self.unqiue_parts has {len(self.unique_parts)} elements. The two must match.")
        parts_dict = self._unique_parts_linkers("part", *parts)
        basic_assemblies = []
        for assembly in self.basic_assemblies:
            parts_linkers = [
                part_linker if isinstance(part_linker, BasicLinker) else parts_dict[_seqrecord_hexdigest(part_linker)]["part"] for part_linker in assembly.parts_linkers
            ]
            basic_assemblies.append(BasicAssembly(assembly.id, *parts_linkers))
        self.__init__(*basic_assemblies)
    
    def _return_clips_data(self):
        """Returns a dictionary of ClipReactions with values describing
        basic_assemblies it uses."""
        clips_dict = OrderedDict(
            **{clip_reaction: [] for assembly in self.basic_assemblies for clip_reaction in assembly.clip_reactions})
        for assembly in self.basic_assemblies:
            for clip_reaction in assembly.clip_reactions:
                clips_dict[clip_reaction].append(assembly)
        return clips_dict

    def _unique_parts_linkers(self, object_key: str, *parts_linkers):
        """Returns a dictionary of unique objects in *parts_linkers. Includes
        an empty list for each item to populate with clip_reactions used by
        each unique part/linker.

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
        """If multiple elements of self.basic_assemblies have same "id"
        attribute, raises a BuildException."""
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

    def export_csv(self):
        """Writes information about each clip and assembly to
        two dependent CSV files in the same folder the command
        is executed"""
        joint_assembly_clips = {}
        for assembly in self.basic_assemblies:
            joint_assembly_clips[assembly.id] = []
        clips = list(self.clips_data.items())
        for index, clip in enumerate(clips):
            for associated_assembly in clip[1]:
                joint_assembly_clips[associated_assembly.id].append(index+1)
        # Writing CSV files
        with open('Clips.csv','w',newline='') as f:
            fieldnames = ['index',
                          'prefix_id',
                          'part_id',
                          'part_name',
                          'suffix_id',
                          'total_assemblies',
                          'assembly_index']
            thewriter = csv.DictWriter(f,fieldnames=fieldnames)
            thewriter.writeheader()
            for item in self._csv_clip_map(joint_assembly_clips):
                thewriter.writerow(item)
        with open('Assemblies.csv','w',newline='') as f:
            fieldnames = ['index','assembly_id','clip_reaction_ids']
            thewriter = csv.DictWriter(f,fieldnames=fieldnames)
            thewriter.writeheader()
            for item in self._csv_assembly_map(joint_assembly_clips):
                thewriter.writerow(item)

    def _csv_clip_map(self,joint_assembly_clips):
        """Returns a List of dictionaries parsable by csv.DictWriter
        of the following format:
        [{'index': 65, 'prefix_id': 'LMP-P', 'part_id': 'B-P63',
          'part_name': 'B-P63_Terminator3_SaITTC_RiboC',
          'suffix_id': 'UTR1-S', 'total_assemblies': 1,
          'assembly_index': [63]},...]

        Args:
            joint_assembly_clips(Dict) -- assembly.id:array of clip indices
            Example -- {'B-P1': [1, 2, 3], 'B-P2': [1, 3, 4],...}
        """
        temp = list(joint_assembly_clips.items())
        Clip_csv_object = []
        for clip in self.clips_data:
            assemblies = []
            for basicassembly in self.clips_data[clip]:
                res = [idx for idx, key in enumerate(temp) if key[0] == basicassembly.id]
                assemblies.append(res[0]+1)
            Clip_csv_object.append({
                'index':list(self.clips_data).index(clip)+1,
                'prefix_id':clip._prefix.prefix_id,
                'part_id':clip._part.id,
                'part_name':clip._part.name,
                'suffix_id':clip._suffix.suffix_id,
                'total_assemblies':len(assemblies),
                'assembly_index':assemblies,
            })
        return Clip_csv_object 

    def _csv_assembly_map(self,joint_assembly_clips):
        """Returns a List of dictionaries parsable by csv.DictWriter
        of the following format:
        [{'index': 1, 'assembly_id': 'B-P1', 'clip_reaction_ids': [1, 2, 3]}
         ,...]

        Args:
            joint_assembly_clips(Dict) -- assembly.id:array of clip indices
            Example -- {'B-P1': [1, 2, 3], 'B-P2': [1, 3, 4],...}
        """
        Assembly_csv_object = []
        for item in joint_assembly_clips:
            Assembly_csv_object.append({
                'index':list(map(lambda assembly: assembly.id,self.basic_assemblies)).index(item)+1,
                'assembly_id':item,
                'clip_reaction_ids':joint_assembly_clips[item],
            })
        print(Assembly_csv_object)
        return Assembly_csv_object

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


class BuildException(Exception):
    pass


def _seqrecord_hexdigest(seqrecord_obj):
        """Returns an MD5 hash of a Bio.SeqRecord.SeqRecord-like object, using
        relevant attributes."""
        seqrec_hash = hashlib.md5(str(seqrecord_obj.seq).encode("UTF-8"))
        bytes_objs = [getattr(seqrecord_obj, attribute).encode("UTF-8") for attribute in [
            "id",
            "name",
            "description"
        ]]
        for element in bytes_objs:
            seqrec_hash.update(element)
        return seqrec_hash.hexdigest()