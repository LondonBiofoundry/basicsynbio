# basicsynbio

An open-source Python API to facilitate BASIC DNA Assembly workflows

[![Build Status][travis_badge]][travis_url]
[![Coverage Status][coverage_badge]][coverage_url]
[![pypi badge][pypi_badge]][pypi_url]
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

## Installation

```sh
pip install basicsynbio
```

## Usage

Import a BASIC part, specified in a genbank file:

```python
import basicsynbio as bsb

basic_part = bsb.import_part("basic_part.gb", "genbank")
```

Or mutliple basic_parts listed seperately in the same file:

```python
basic_parts = bsb.import_parts("basic_parts.gb", "genbank")
```

Make an assembly using several BASIC parts and any supplied [Biolegio linkers](https://www.biolegio.com/products-services/basic/), then return an annotated genbank file:

```python
assembly = bsb.BasicAssembly(
    bsb.biolegio_dict["LMP"],
    basic_parts[0],
    bsb.biolegio_dict["LMS"],
    basic_parts[1]
)
assembly.return_file("my_basic_assembly.gb")
```

## Contributing

1. Fork it <https://github.com/hainesm6/basicsynbio/fork>
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -m 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

## Meta

This project is licensed under the MIT License - see the ``LICENSE`` file for details

Originally created by `Matthew Haines - hainesm6@gmail.com`

[pypi_badge]: https://img.shields.io/pypi/v/basicsynbio.svg
[pypi_url]: https://pypi.python.org/pypi/basicsynbio
[travis_badge]: https://travis-ci.org/hainesm6/basicsynbio.svg
[travis_url]: https://travis-ci.org/hainesm6/basicsynbio
[coverage_badge]:
[coverage_url]:
