Part and linker collections
===========================

Introduction
------------

As described and illustrated in the usage section (:ref:`accessing-part-collections`), basicsynbio contains a collection of Parts and Linkers accessible 
from within the API. In contrast to previous approaches, accessing these collections doesn't require 3rd-party database access or even an internet connection.
Below is a list of selected collections to browse and instructions for contributing part collections to the API.

.. _browse-collections:

Browse selected collections
---------------------------

.. note::

    Please print the following collections for information on them or see the web-app if available: ``BASIC_BIOLEGIO_LINKERS["v0.1"]``, ``BASIC_SEVA_PARTS["v0.1"]``, ``BASIC_PROMOTER_PARTS["v0.1"]``

.. todo::

    Add links to printed versions of the collections and images of plasmid maps.

Contributing part collections to basicsynbio
--------------------------------------------

The following steps are required to add new part collections or versions of part collections to the API.
*Please feel free to contact* `hainesm6`_ *directly for help with this*:

.. _hainesm6: mailto:hainesm6@gmail.com

#. Generate a genbank file containing all parts in the collection as separate elements.
#. Via a pull-request (refer to :doc:`contributing`), complete the following:

   * Add the genbank file to the :py:mod:`babasicsynbio.parts_linkers` sub-package.
   * Using :py:mod:`importlib.resources` and :py:func:`basicsynbio.parts_linkers.main.make_collection`, make or update the collection in :py:mod:`basicsynbio.parts_linkers.basic_parts`.
   * Where required, update the parts_linkers sub-package \_\_init\_\_.py file to import the collection.

