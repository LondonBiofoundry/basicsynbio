# basicsynbio

An open-source Python API to facilitate [BASIC DNA Assembly](https://www.basic-assembly.org/about) workflows

[![Build Status][travis_badge]][travis_url]
[![Coverage Status][coverage_badge]][coverage_url]
[![pypi badge][pypi_badge]][pypi_url]
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

## Installation

*Available on PyPi soon. For now clone repo and install in [development mode][development_url] or build and install wheels (refer to [PyPA](https://www.pypa.io/en/latest/)).*

## Usage

*The basicsynbio API extends the [Bioython library](https://biopython.org/). Extensive knowledge of Biopython is not required but basic knowledge of key objects would aid users.*

To create genbank or FASTA files for BASIC DNA assemblies:

1. [Import BASIC parts](#1-import-basic-parts)
2. [Create assemblies from BASIC parts and linkers](#2-create-assemblies-from-basic-parts-and-linkers)
3. [Return new parts and hierarchical assemblies](#3-return-new-parts-and-hierarchical-assemblies)
4. [Exporting sequences](#4-exporting-sequences)

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

All BASIC parts require flanking iP and iS sequences. To add these use the optional `add_i_seqs` argument, available for all the above functions e.g.

```python
basic_part = bsb.seqrec2part(SeqRecord, add_i_seqs=True)
```

### 2. Create assemblies from BASIC parts and linkers

Make an assembly using several BASIC parts and any [Biolegio linkers](https://www.biolegio.com/products-services/basic/), then return an annotated genbank file:

```python
assembly = bsb.BasicAssembly(
    bsb.biolegio_dict["LMP"],
    basic_parts[0],
    bsb.biolegio_dict["LMS"],
    basic_parts[1]
)
```

### 3. Return new parts and hierarchical assemblies

The `return_part()` method returns a BasicPart object from the assembly, enabling downstream BasicAssemblies (*LMP and LMS flank new part sequences*):

```python
new_part = assembly.return_part(id="new_part")
hierarchical_assembly = bsb.BasicAssembly(
    new_part,
    ...
)
```

### 4. Exporting sequences

Assemblies can be exported to a file directly:

```python
assembly.return_file(handle="hello_world_assembly.gb")
```

or if necessary a Biopython SeqRecord created:

```python
seqrec = assembly.return_seqrec(id="assembly")
```

The `export_to_file()` function exports several objects to the same file:

```python
bsb.export_to_file(
    basic_parts,
    handle="another_basic_parts_file.gb",
)
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

