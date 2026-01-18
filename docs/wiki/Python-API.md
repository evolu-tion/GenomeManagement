# Python API Reference

You can use the internal classes of GenomeManagement for custom scripts.

## `FastaManager`

Handles FASTA file operations and sequence manipulation.

```python
from genomemanagement.seq_manage import FastaManager

# Initialize
fasta = FastaManager("genome.fa", show_genome_stat=False)

# Get Sequence
seq = fasta.get_sequence("Chr1", 100, 200, "+")

# Check Chromosome
exists = fasta.check_chromosome("Chr1")

# Utilities
rev_comp = fasta.complementary("ATGC")
gc_content = fasta.get_gc_content("ATGC")
```

## `GffManager`

Parses and queries GFF3 annotation files.

```python
from genomemanagement.seq_manage import GffManager

# Initialize
gff = GffManager("annotations.gff")

# Check Gene
exists = gff.check_gene("GeneID")

# Get Gene Structure
structure = gff.get_table_data_of_gene("GeneID")
```

## `GenomeManager`

Combines FASTA and GFF functionality for advanced operations. Inherits from both `FastaManager` and `GffManager`.

```python
from genomemanagement.seq_manage import GenomeManager

# Initialize
genome = GenomeManager("genome.fa", "annotations.gff")

# Advanced Retrieval
# (See source for complex promoter retrieval methods)
```
