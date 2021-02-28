import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="basicsynbio",
    version="0.4.0",
    url="https://github.com/LondonBiofoundry/basicsynbio",
    license="BSD-3-Clause License",
    author="LondonBiofoundry",
    author_email="hainesm6@gmail.com",
    description="An open-source Python package to facilitate BASIC DNA Assembly workflows",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(
        exclude=("tests")
    ),
    include_package_data=True,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
    install_requires=[
        "biopython>=1.78",
        "pandas",
        "platemap",
        "python-Levenshtein",
        "reportlab",
        "sbol2",
    ],
    project_urls={
        "Documentation": "https://londonbiofoundry.github.io/basicsynbio/index.html",
        "Source": "https://github.com/LondonBiofoundry/basicsynbio",
    },
)
