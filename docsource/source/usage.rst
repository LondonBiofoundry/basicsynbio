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

.. _accessing-part-collections:

1a. Accessing BASIC part collections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* The basicsynbio package contains collections of parts and linkers compatible with BASIC DNA assembly. 
* Each collection can have multiple versions e.g. ``"v0.1"``, with changes only in minor numbers indicating backwards compatibility e.g. ``"v0.1"`` would be compatible with ``"v0.2"``.
* Within each version of a collection are individual part or linker objects. For instance, to access the BASIC backbone with ampicilin resistance and a pUC ori (equivalent to `SEVA`_ 18), input the following:

.. _SEVA:  <http://seva-plasmids.com/>

.. code-block:: python

    import basicsynbio as bsb
    basic_seva18 = bsb.BASIC_SEVA_PARTS["v0.1"]["18"]

* A list of selected part and linker collections is given in :ref:`browse-collections`.
* The contents of each collection can also be displayed using the print function e.g.

.. code-block:: python

    print(bsb.BASIC_SEVA_PARTS["v0.1"])

1b. Import parts from external sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compatible BASIC parts can be imported from multiple sources.
To import one BASIC part from a local file (e.g. genbank file) run:

.. code-block:: python

    import basicsynbio as bsb

    basic_part = bsb.import_part("basic_part.gb", "genbank")

Otherwise mutliple BASIC parts from the same gb file are imported as follows:

.. code-block:: python

    basic_parts = bsb.import_parts("basic_parts.gb", "genbank")

Convert `Biopython SeqRecords`_ or similar objects into BASIC parts.
This is useful for accessing sequences via `NCBI Entrez`_ directly. 
The following leverages BioPython and Entrez to generate a new BASIC part encoding sfGFP from `KJ541673.2`_:

.. _Biopython SeqRecords: https://biopython.org/wiki/SeqRecord
.. _NCBI Entrez: https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html
.. _KJ541673.2: https://www.ncbi.nlm.nih.gov/nuccore/KJ541673.2

.. code-block:: python

    from Bio import Entrez, SeqIO
    import basicsynbio as bsb
    from basicsynbio.utils import feature_from_qualifier
    Entrez.email = "hainesm6@gmail.com"
    with Entrez.efetch(db="nucleotide", id="KJ541673.2", rettype="gb", retmode="text") as handle:
        kj541673 = SeqIO.read(handle, "genbank")
        sfgfp_feature = feature_from_qualifier(kj541673, "gene", ["sfGFP"])
        sfgfp = kj541673[sfgfp_feature.location.start:sfgfp_feature.location.end]
    sfgfp_part = bsb.seqrec2part(sfgfp, add_i_seqs=True)

For parts specified in `SBOL <https://sbolstandard.org/>`_ the following imports them as a generator object:

.. code-block:: python

    basic_parts = bsb.import_sbol_part("basic_parts.rdf")

All BasicPart objects require flanking *i*\ P and *i*\ S sequences. To add these
when creating your object, use the optional ``add_i_seqs`` argument,
available for all the above functions e.g.

.. code-block:: python

    basic_part = bsb.seqrec2part(SeqRecord, add_i_seqs=True)

2. Create the assemblies
~~~~~~~~~~~~~~~~~~~~~~~~

Create a ``BasicAssembly`` object from your imported BASIC parts using any
`Biolegio Linkers`_ contained within the ``BIOLEGIO_LINKERS`` collection:

.. _Biolegio Linkers: https://www.biolegio.com/products-services/basic/ 

.. code-block:: python
    
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
subsequent hierarchical assembly. Use the ``return_part()`` method on a BasicAssembly object to simulate this behaviour:

.. code-block:: python

    new_part = assembly.return_part(name="new part from assembly")
    hierarchical_assembly = bsb.BasicAssembly(
        new_part,
        ...
    )

3. Create the build
~~~~~~~~~~~~~~~~~~~

More often than not, a collection of BASIC assemblies are constructed in parallel. 
To aid this process users should create a ``BasicBuild`` object using multiple
BasicAssembly objects:

.. literalinclude:: /literal_includes/make_literals.py
    :pyobject: build_json
    :start-after: def build_json():
    :end-before: return build
    :dedent: 4

The ``build`` instance contains data describing the unique BasicParts, BasicLinkers and ClipReactions objects
associated with this build. These objects together provide a description of the materials and steps required
to construct your assemblies.

4. Export your data
~~~~~~~~~~~~~~~~~~~

BasicBuild objects can be serialised as JSON or
exported as two zipped csv files describing build clip reactions and assemblies.

The ``export_csvs()`` method generates a zip file, containing :doc:`clips_csv` and :doc:`assemblies_csv`:

.. code-block::
    
    build.export_csvs("build_csvs.zip")

The ``bsb.pdf_instructions()`` function creates pdf instructions for manual assembly of the build in the lab. An example can be seen here_

.. _here: https://github.com/LondonBiofoundry/basicsynbio/blob/master/docsource/source/literal_includes/export_pdf_example.pdf

.. code-block::
    
    bsb.export_csvs(build)

The ``bsb.export_echo_assembly(build)`` function creates echo lab automation instructions for the clips to assembly step of the build. Example outputs, :doc:`echo_clips_1_csv` and :doc:`echo_water_buffer_1_csv`

.. code-block::
    
    bsb.export_echo_assembly(build)

To serialise the build, the `json API`_ can be used, in the following case yielding (:doc:`build_json`):

.. _json API: https://docs.python.org/3/library/json.html

.. literalinclude:: /literal_includes/make_literals.py
    :pyobject: export_json
    :start-after: def export_json(build):
    :dedent: 4

Depending on the file format, the resulting output 
contains data on the unique BasicParts (``unique_parts``), BasicLinkers (``unique_linkers``)
and ClipReactions (``clips_data``) objects required to build the assemblies (``assembly_data``).
This data can either be analysed directly, informing manual workflows or further processed to 
generate arguments for liquid-handling systems.

In addition to exporting build data, **it is recommended to export
annotated BasicAssembly objects and the unique BasicParts** associated with the build. The later is important
for completly decoding serialised BasicBuild objects, described in the next section.
Notably, any collection of BasicPart [#f1]_ or BasicAssembly
objects can be exported using the formats supported by `BioPython`_, with the default being genbank:

.. _BioPython: https://biopython.org/wiki/SeqIO

.. code-block:: python

    unique_parts = build.unique_parts
    bsb.export_sequences_to_file(unique_parts, "the_parts_i_need.gb")
    bsb.export_sequences_to_file(promoter_assemblies, "promoter_assemblies.gb")

Importing from build.json
-------------------------

It is possible to decode build.json objects, restoring the BasicBuild object.
Users have two options:

#. The first method uses only the build.json file and results in correct sequences, although, with a loss of metainformation e.g. annotations, features etc.
#. The second method extends the first, updating the decoded BasicBuild object using the original BasicParts with the correct annotations. 

To partially decode a build.json file:

.. code-block:: python

    import basicsynbio as bsb
    import json

    with open("build.json") as build_json:
        partially_decoded_build = json.load(build_json, cls=bsb.BuildDecoder)

To completely decode this file:

.. code-block:: python

    original_parts = bsb.import_parts("the_parts_i_need.gb", "genbank")
    decoded_build.update_parts(*original_parts)

.. rubric:: Footnotes

.. [#f1] This also applies to any Biopython SeqRecord-like object.
