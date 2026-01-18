# Command Line Tools

GenomeManagement provides several command-line tools for common tasks.

## 1. promoter-retrieve

Retrieves promoter sequences from a genome based on gene annotations.

**Usage:**

```bash
promoter-retrieve \
    --genome <genome.fa> \
    --gff <annotations.gff> \
    --output <output.fa> \
    --output_format <fasta/gff> \
    --type <TSS/TLS> \
    --upstream <bp> \
    --downstream <bp>
```

**Options:**

- `-g, --genome`: Input genome FASTA file.
- `-f, --gff`: Input GFF3 annotation file.
- `-o, --output`: Output file path.
- `-p, --output_format`: Output format (`fasta` or `gff`).
- `-t, --type`: Reference point (`TSS` = Transcription Start Site, `TLS` = Translation Start Site).
- `-u, --upstream`: Base pairs upstream of reference.
- `-d, --downstream`: Base pairs downstream of reference.
- `-a, --all_gene`: Retrieve for all genes (`Yes`/`No`). Default: `Yes`.
- `-l, --selected_gene_list`: Text file with list of gene IDs (required if `--all_gene No`).

## 2. proteins-retrieve

Retrieves protein sequences for a specific list of protein IDs.

**Usage:**

```bash
proteins-retrieve \
    --input <proteins.fa> \
    --list_of_interest <id_list.txt> \
    --output <output.fa>
```

## 3. get-genome-statistic

Calculates and displays statistics for a genome assembly (N50, GC content, etc.).

**Usage:**

```bash
get-genome-statistic --genome <genome.fa>
```

## 4. plantpan2 & plantpan3

Automated tools for retrieving promoter analysis data from PlantPAN databases.

**Usage:**

```bash
plantpan2
# OR
plantpan3
```

_Note: These scripts may require specific input file locations or interactive configuration as per the original script design._

## 5. meme-xml2bed

Converts MEME motif search results (XML) to BED format.

**Usage:**

```bash
meme-xml2bed \
    --input <meme.xml> \
    --output <meme.bed> \
    --pvalue <1e-2>
```
