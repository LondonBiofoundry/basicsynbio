# Usage

*The basicsynbio API extends the [Bioython library](https://biopython.org/). Extensive knowledge of Biopython is not required but basic knowledge of key objects would aid users.*

Assemble BASIC DNA assembly constructs and export the assembled sequences to FASTA or GenBank files:

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
bsb.export_sequences_to_file([assembly1, part1], "file.gb")