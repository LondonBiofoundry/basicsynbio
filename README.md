# basicsynbio

An open-source Python package to facilitate [BASIC DNA Assembly](https://www.basic-assembly.org/about) workflows

[![Build Status][travis_badge]][travis_url]
[![Coverage Status][coverage_badge]][coverage_url]
[![pypi badge][pypi_badge]][pypi_url]
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

## Installation

*Available on PyPi soon. For now clone repo and install in [development mode][development_url] or build and install wheels (refer to [PyPA](https://www.pypa.io/en/latest/)).*

## Usage

*The basicsynbio API extends the [Bioython library](https://biopython.org/). Extensive knowledge of Biopython is not required but basic knowledge of key objects would aid users.*

Assemble BASIC DNA assembly constructs and export the assembled sequences to FASTA or GenBank files:

- [basicsynbio](#basicsynbio)
  - [Installation](#installation)
  - [Usage](#usage)
    - [1. Import BASIC parts](#1-import-basic-parts)
    - [2. Create assemblies from BASIC parts using linkers](#2-create-assemblies-from-basic-parts-using-linkers)
    - [3. Create hierarchical assemblies](#3-create-hierarchical-assemblies)
    - [4. Export sequences](#4-export-sequences)
  - [Contributing](#contributing)
  - [Meta](#meta)

### 1. Import BASIC parts

Import a BASIC part from a file (e.g. genbank file):

```python
import basicsynbio as bsb

basic_part = bsb.import_part("basic_part.gb", "genbank")
```

Or mutliple basic_parts listed seperately in the same file:

```python
basic_parts = bsb.import_parts("basic_parts.gb", "genbank")
```

Alternatively, convert a Biopython SeqRecord object into a BASIC part:

```python
basic_part = bsb.seqrec2part(SeqRecord)
```

All BasicPart objects require flanking iP and iS sequences. To add these when creating your object, use the optional `add_i_seqs` argument, available for all the above functions e.g.

```python
basic_part = bsb.seqrec2part(SeqRecord, add_i_seqs=True)
```

### 2. Create assemblies from BASIC parts using linkers

Create a BasicAssembly object from your imported BASIC parts using any [Biolegio linkers](https://www.biolegio.com/products-services/basic/) contained within the `biolegio_dict` dictionary:

```python
assembly = bsb.BasicAssembly(
    bsb.biolegio_dict["LMP"],
    basic_parts[0],
    bsb.biolegio_dict["LMS"],
    basic_parts[1]
)
```

### 3. Create hierarchical assemblies

The `return_part()` method returns a BasicPart from the assembly, enabling hierarchical assemblies. *LMP and LMS linkers must flank the new part sequence within the initial assembly*:

```python
new_part = assembly.return_part(id="new_part")
hierarchical_assembly = bsb.BasicAssembly(
    new_part,
    ...
)
```

### 4. Export sequences

Export `BasicPart`, `BasicAssembly` and `Bio.SeqRecord` objects to a file. One or more can be given:

```python
bsb.export_to_file([assembly1, part1], "file.gb")
```

## Contributing

1. Fork it <https://github.com/hainesm6/basicsynbio/fork>
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -m 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

## Meta

This project is licensed under the MIT License - see the ``LICENSE`` file for details

[pypi_badge]: https://img.shields.io/pypi/v/basicsynbio.svg
[pypi_url]: https://pypi.python.org/pypi/basicsynbio
[travis_badge]: https://travis-ci.org/hainesm6/basicsynbio.svg
[travis_url]: https://travis-ci.org/hainesm6/basicsynbio
[coverage_badge]: https://coveralls.io/repos/github/hainesm6/basicsynbio/badge.svg?branch=master
[coverage_url]: https://coveralls.io/github/hainesm6/basicsynbio?branch=master
[development_url]: https://packaging.python.org/guides/distributing-packages-using-setuptools/#working-in-development-mode
