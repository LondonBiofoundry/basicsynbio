# json representation of BASIC DNA assembly builds

This document tries to come with a json representation of BASIC DNA assembly builds.

## Requirements for json representation

1. Object easily converted into build instructions for liquid-handling robotics.
2. *Object facilitates validation of constructs.*
3. Capable of decoding object back into a BasicBuild object.

## Draft data structure

```json
build = {
    "unique_parts": {
        "hash((part_1.id, part_1.seq))": {
            "id": "part_1.id",
            "sequence": "part_1.seq",
            "clips_indexes": []
        }, ...
    },
    "unique_linkers": {
        "hash((linker_1.id, linker_1.seq))": {
            "id": "linker_1.id",
            "linker_class": "linker_1 class",
            "sequence": "linker_1.seq",
            "prefix_id": "linker_1 prefix",
            "suffix_id": "linker_1 suffix",
            "clips_indexes": []
        }, ...
    },
    "clips_data": [
        {
            "prefix_linker_key": "hash of corresponding linker from unique linkers",
            "part_key": "hash of correspond part from unique parts",
            "suffix_linker_key": "hash of corresponding linker from unique linkers",
            "basic_assemblies": [
                "index corresponding to 1st basic assembly using this clip reaction", ...
            ]
        }, ...
    ],
    "basic_assemblies": [
        {
            "id": "basic_assembly_1.id",
            "clips_indexes": []
        }, ...
    ],
}
```

## How build json object facilitates instructions for building assemblies?

Specific build instructions would centre around setting up clip reactions and completing assemblies. Arguments for purification and transformation steps are easily inferred from arguments used to setup these two critical processes. 

To setup clip reactions:
- `"clips_data"` contains information on what prefix, suffix and part to transfer including id and sequence. This information could be used to query databases describing collections of parts/linkers in lab freezers, enabling substrates to be parsed into liquid-handling robotic jobs.
- Parts sequences in combination with concentration values can be used to ensure the correct volume of part is transfered during clip reaction setups.
- To calculate the number of times a specific clip reaction must be repeated in a build, it provides the `"basic_assemblies"` array which lists which assemblies use this specific clip reaction. The number of times the i<sup>th</sup> clip reaction must be repeated is calculated by:

```python
number = len(build.clips_data[i].basic_assemblies)/assemblies_per_clip
```

where `assemblies_per_clip` will vary depending on the liquid-handling platform used e.g. Opentrons OT-2, Lacyte Echo etc.

To complete assemblies:
- Each object in `build.basic_assemblies` contains the indexes of required clip reactions as an array (`clips_indexes`). These can be used to identify which specific clip reactions to transfer when assembling the construct.

## How to validate constructs from a build?

Using the basicsynbio package, export BasicAssembly objects within the `BasicBuild.basic_assemblies` attribute to a GenBank file using the `export_sequences_to_file` function:

```python
import basicsynbio as bsb
my_build = bsb.BasicBuild(*my_assemblies)
bsb.export_sequences_to_file(my_build.basic_assemblies, handle=file_handle)
```

The id attribute of a `build.basic_assemblies` element will match the [VERSION field](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#VersionB) of the associated genbank entry in the resulting file. A  Basic Assembly construct has been successfully assembled if it's sequence matches that prediceted by the genbank file. Techniques such as Sanger Sequencing and diagnostic digests can be used to facilitate this (Ref).

*Can then specifically mention these approaches. Specifically using sequencing primers that anneal to the T0 & T1 for validating inserts and BsaI diagnostic digests for validating insert/backbone*.

## How the build json object is encoded/decoded?

BasicBuild objects can be encoded as follows:

```python
import basicsynbio as bsb
my_build = bsb.BasicBuild(*my_assemblies)
json.dumps(my_build, cls=bsb.BuildEncoder)
```

To prevent loss of information (discussed later), it is recommended to retain BasicBuild.unique_parts object data, as a separate genbank file (below) or via pickling:

```python
bsb.export_sequences_to_file(my_build.unique_parts, file_handle)
```

There are two options for decoding build.json objects:

1. The first method uses only the build.json object and results in correct sequences, although, with a loss of metainformation e.g. annotations, features etc.
2. The second method uses the build.json object in addition to a collection of BasicPart objects equivalent to that stored in the BasicBuild.unqiue_parts attribute at the time of serialisation. This completely recapitulates the BasicBuild object at the time of serialisation.

An example of the 1st method is given by:

```python
import basicsynbio as bsb
partially_decoded_build = json.loads(build_json, object_hook=bsb.build_object_hook)
```

An example of the 2nd method is given by:

```python
import basicsynbio as bsb
build_unique_parts = bsb.import_parts(build_unqiue_parts_handle, "genbank")
decoded_build = bsb.decode_build(build_json, *build_unique_parts)
```

