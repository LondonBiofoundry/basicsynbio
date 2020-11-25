Usage
=====

The basicsynbio API extends the `Biopython library <https://biopython.org/>`_.
Extensive knowledge of Biopython is
not required but basic knowledge of key objects would aid users.

basicsynbio workflow
--------------------

The core basicsynbio workflow has the following steps associated with
it:

#. Get the parts (either from collections or externally)
#. Create the assemblies
#. Create the build
#. Export your data

1a. Accessing BASIC part collections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The basicsynbio package contains collections of parts and linkers compatible with BASIC DNA assembly. 
* Each collection behaves like a `dictionary`_ object.
* For instance, to access the BASIC backbone with ampicilin resistance and a pUC ori (equivalent to `SEVA`_ 18), input the following:

.. _dictionary: <https://docs.python.org/3/tutorial/datastructures.html#dictionaries>
.. _SEVA:  <http://seva-plasmids.com/>

.. code:: python

    import basicsynbio as bsb
    basic_seva18 = bsb.BASIC_SEVA_PARTS["v0.1]["18"]

* A list of all part and linker collections is given :doc:`collections`.
* The contents of each collection can be displayed using the print function e.g.

.. code:: python

    print(bsb.BASIC_SEVA_PARTS["v0.1])

1b. Import parts from external sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compatible BASIC parts can be imported from multiple sources.
To import one BASIC part from a local file (e.g. genbank file) run:

.. code:: python

    import basicsynbio as bsb

    basic_part = bsb.import_part("basic_part.gb", "genbank")

Otherwise mutliple BASIC parts from the same gb file are imported as follows:

.. code:: python

    basic_parts = bsb.import_parts("basic_parts.gb", "genbank")

Convert `Biopython SeqRecords`_ or similar objects into BASIC parts:

.. _Biopython SeqRecords: https://biopython.org/wiki/SeqRecord

.. code:: python

    basic_part = bsb.seqrec2part(SeqRecord)

For a part given in `SBOL <https://sbolstandard.org/>`_ use:

.. code:: python

    basic_part = bsb.import_sbol_part("basic_part.rdf")

Finally you can also import one or more parts from a JBEI-ICE instance, e.g. the `public-registry`_:

.. _public-registry: https://public-registry.jbei.org/

.. code:: python

    ice_nums = (string(int) for int in range(17297, 17339))
    basic_parts = bsb.import_ice_parts(ice_user_config, *ice_nums)  


All BasicPart objects require flanking *i*\ P and *i*\ S sequences. To add these
when creating your object, use the optional ``add_i_seqs`` argument,
available for all the above functions e.g.

.. code:: python

    basic_part = bsb.seqrec2part(SeqRecord, add_i_seqs=True)

2. Create the assemblies
~~~~~~~~~~~~~~~~~~~~~~~~

Create a ``BasicAssembly`` object from your imported BASIC parts using any
`Biolegio Linkers`_ contained within the ``BIOLEGIO_LINKERS`` collection:

.. _Biolegio Linkers: https://www.biolegio.com/products-services/basic/ 

.. code:: python
    
    import basicsynbio as bsb
    my_basic_part = bsb.import_part("my_basic_part.gb", "genbank")
    assembly = bsb.BasicAssembly(
        "my_first_basic_assembly",
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
        my_basic_part,
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"],
        bsb.BASIC_SEVA_PARTS["v0.1"]["18"]
    )

This creates a BasicAssembly object where ``my_basic_part`` has been cloned
into the BASIC_SEVA_18 backbone.

A desirable feature of BASIC DNA Assembly is its single-tier format (:doc:`introduction`).
This ensures any assembly flanked by LMP and LMS linkers can be used in a 
subsequent hierarchical assembly:

.. code:: python

    new_part = assembly.return_part(id="new_part")
    hierarchical_assembly = bsb.BasicAssembly(
        new_part,
        ...
    )

3. Create the build
~~~~~~~~~~~~~~~~~~~

More often than not, a collection of BASIC assemblies are constructed in parallel. 
To aid this process users should create a ``BasicBuild`` object using multiple
BasicAssembly objects:

.. code:: python

    import basicsynbio as bsb

    promoter_assemblies = (bsb.BasicAssembly(
        f"promoter_construct_{ind}",
        bsb.BASIC_SEVA_PARTS["v0.1"]["26"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMP"],
        promoter,
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["UTR1-RBS2"],
        bsb.BASIC_CDS_PARTS["v0.1"]["sfGFP"],
        bsb.BASIC_BIOLEGIO_LINKERS["v0.1"]["LMS"]
        ) for ind, promoter in enumerate(bsb.BASIC_PROMOTER_PARTS["v0.1"].values()))
    build = bsb.BasicBuild(*promoter_assemblies)

The ``build`` instance contains data describing the unique BasicParts, BasicLinkers and ClipReactions objects
associated with this build. These objects together provide a description of the materials and steps required
to construct your assemblies.

4. Export your data
~~~~~~~~~~~~~~~~~~~

BasicBuild objects can be serialised using the `json API`_, part of the standard library:

.. _json API: https://docs.python.org/3/library/json.html

.. code:: python
    
    import json

    with open("my_build.json", "w") as json_file:
        json.dump(build, json_file, cls=bsb.BuildEncoder, indent=4)

Like the associated build object, the resulting output (:doc:`build_example`)
contains data on the unique BasicParts (``unique_parts``), BasicLinkers (``unique_linkers``)
and ClipReactions (``clips_data``) objects required to build the assemblies (``assembly_data``).
This data can either be analysed directly or further processed to 
generate assemblies manually or via any liquid-handling robotic platform.

In addition to exporting build data as a json file, **it is recommended to export
annotated BasicAssembly objects and the unique BasicParts** associated with the build.
Notably, any collection of BasicPart [#f1]_ or BasicAssembly
objects can be exported using the formats supported by `BioPython`_:

.. _BioPython: https://biopython.org/wiki/SeqIO

.. code:: python

    unique_parts = (part_dict["part"] for part_dict in build.unique_parts.values())
    bsb.export_sequences_to_file(unique_parts, "the_parts_i_need.gb")
    bsb.export_sequences_to_file(cds_assemblies, "cds_assemblies.gb")

Importing from build.json
-------------------------

It is possible to decode build.json objects, restoring the BasicBuild object.
Users have two options:

#. The first method uses only the build.json file and results in correct sequences, although, with a loss of metainformation e.g. annotations, features etc.
#. The second method extends the first, updating the decoded BasicBuild object using the original BasicParts with the correct annotations. 

To partially decode a build.json file:

.. code:: python

    import basicsynbio as bsb
    import json

    with open("build.json") as build_json:
        partially_decoded_build = json.load(build_json, cls=bsb.BuildDecoder)

To completely decode this file:

.. code:: python

    original_parts = bsb.import_parts("the_parts_i_need.gb", "genbank")
    decoded_build.update_parts(*original_parts)

.. rubric:: Footnotes

.. [#f1] This also applies to any Biopython SeqRecord-like object.
