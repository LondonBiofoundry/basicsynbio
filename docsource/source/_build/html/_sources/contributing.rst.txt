Contributing
============

Setting up a development environment
------------------------------------

#. Fork the `basicsynbio repository`_.
#. Clone your fork::

    git clone https://github.com/username/basicsynbio.git

#. cd into the directory.
#. Using a virtual environment running Python 3.8 or greater install the dependencies::

    pip install -r requirements.txt

#. Install the package in editable mode::

    pip install -e .

#. Confirm all test pass::

    pytest [--runslow]

.. _basicsynbio repository: https://github.com/LondonBiofoundry/basicsynbio.git

Making changes
--------------

Steps
^^^^^

#. set the upstream repository::

    git remote add upstream https://github.com/LondonBiofoundry/basicsynbio.git

#. Pull the latest changes to origin/master::

    git pull upstream master

#. Create a new branch for your development::

    git checkout -b dev/foobar

#. Make the development then stage and commit::

    git add .
    git commit -m "make this desirable development to basicsynbio"

#. Push the development back to GitHub::

    git push origin dev/foobar

#. Go to GitHub and create a new Pull Request.

Guidelines for changes
----------------------

General Guidelines
^^^^^^^^^^^^^^^^^^

* Code should be formatted using `black`_.
* New code should have associated tests written using `pytest`_.
* Update the documentation if required. The documentation is formatted using `Sphinx`_ and as such all files must be written in reST.
* Update GitHub pages by running ``make github`` in the */docsource* directory. This implements GitHub Pages as described by `Anne Gentle`_.
* New releases are made online (https://github.com/LondonBiofoundry/basicsynbio/releases). Releases are automatically published to PyPi via the `build-publish.yml workflow`_.

.. _Anne Gentle: https://www.docslikecode.com/articles/github-pages-python-sphinx/
.. _build-publish.yml workflow: https://github.com/LondonBiofoundry/basicsynbio/blob/master/.github/workflows/build-publish.yml

Guidelines for docstrings
^^^^^^^^^^^^^^^^^^^^^^^^^

* Format docstrings according to the `Google Style Guide`_. This enables automated documentation of the API using `sphinx.ext.napoleon`_.
* Where the same argument is used by multiple methods/functions, the :py:func:`basicsynbio.decorators.addargs2docs` decorator may be used. Associated instances of :py:class:`basicsynbio.decorators.ArgDescription` are placed in :py:class:`basicsynbio.main.CommonArgDocs`. 


.. _black: https://github.com/psf/black
.. _pytest: https://docs.pytest.org/en/stable/
.. _Google Style Guide: https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings
.. _sphinx.ext.napoleon: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html
.. _Sphinx: https://www.sphinx-doc.org/en/master/usage/quickstart.html

