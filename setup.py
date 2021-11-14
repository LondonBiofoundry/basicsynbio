import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

# __version__ is overwritten by GitHub actions workflow `build-publish.yml`
__version__ = "develop"

setuptools.setup(
    name="basicsynbio",
    version=__version__,
    url="https://github.com/LondonBiofoundry/basicsynbio",
    license="BSD-3-Clause License",
    author="LondonBiofoundry",
    author_email="hainesm6@gmail.com",
    description="An open-source Python package to facilitate BASIC DNA Assembly workflows",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(exclude=("tests")),
    include_package_data=True,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.78",
        "pandas",
        "platemap",
        "primer3-py",
        "python-Levenshtein",
        "reportlab",
        "sbol2",
    ],
    project_urls={
        "Documentation": "https://londonbiofoundry.github.io/basicsynbio/index.html",
        "Source": "https://github.com/LondonBiofoundry/basicsynbio",
    },
)
