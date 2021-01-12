Introduction
============

This section provides some background to basicsynbio.

BASIC DNA Assembly
------------------

BASIC DNA Assembly uses DNA parts flanked by standardised *i*\ P and *i*\ S sequences,
along with linker sequences to assemble multi-part DNA constructs. 
It is a powerful assembly technology with > 90 % accuracy for assembling DNA constructs containing 
up to 14 parts and linkers per assembly round [#Storch2015]_. BASIC also benefits from
a single-tier format, where any BASIC part can be used in any BASIC assembly. This
simplifies the workflow, making BASIC amenable to automation [#Storch2020]_ and aiding hierarchial assemblies.
Potential applications include the assembly of pathways, circuits and the generation of fusion proteins.

Why basicsynbio was needed?
---------------------------

As described in the previous section, BASIC DNA assembly has several desirable features.
However, there were a number of issues that kept cropping up for new and seasoned users.

I can't access Part and Linker sequences...
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It was often the case that data describing core Parts and Linkers was not in a central location.
This subsequently made it difficult for users to *in silico* assemble their constructs as 
DNA sequences had to first be acquire from multiple sources/users. 

I can't figure out the sequence of my construct/s...
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once a user had a acquired the parts they needed for their assembly, 
there was the challenge of figuring out what the resulting sequence of the construct was.
This step is critical in determining whether any assembly had been successful. To the best of our
knowledge, no software was previously available to analyse BASIC parts and linkers
and stitch them together. A user would often have to do all of this manually.

What clip reactions do I need to setup and then how do I mix them?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This was a common question amongst users. Particularly when assembling many constructs.
Although, some work was previously done to address this issue [#Storch2020]_, the
result was specific to certain liquid-handling platforms rather than being generic.


My assembly hasn't worked...
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

New users coming to BASIC DNA assembly were often not familiar with the various requirements
of this technology. For instance:

* They wouldn't flank their parts with *i*\ P and *i*\ S sequences.
* Their part would contain an internal BsaI restriction site.
* Parts would be too small in length, leading to loss during purification.
* Parts and linkers wouldn't alternate.
* They would use the same Linker or Linker-half multiple times in the same assembly.

It would then take some time to identify the fault which often involved the help of experts.

Summary
-------

While BASIC DNA assembly has several desirable features, without good software to supplement it, it has a number of shortcomings
making it difficult for new users to adopt and for existing users to implement.
We built the basicsynbio package with this in mind.

.. [#Storch2015] Storch, M., Casini, A., Mackrow, B., Fleming, T., Trewhitt, H., Ellis, T., & Baldwin, G. S. (2015). BASIC: A New Biopart Assembly Standard for Idempotent Cloning Provides Accurate, Single-Tier DNA Assembly for Synthetic Biology. ACS Synthetic Biology, 4(7), 781–787. https://doi.org/10.1021/sb500356d
.. [#Storch2020] Storch, M., Haines, M. C., & Baldwin, G. S. (2020). DNA-BOT: a low-cost, automated DNA assembly platform for synthetic biology. Synthetic Biology, 5(1). https://doi.org/10.1093/synbio/ysaa010