# Genome Management Package

This tool is copyright 2016 by Nattawet Sriwichai,
King Mongkut's University Technology Thonburi (Bioinformatics and Systems Biology Program), Thailand.
All rights reserved. See the licence text below.

# Installation

## Install from PyPI

    pip install genomemanagement

## Install from GitHub

    pip install git+https://github.com/evolu-tion/GenomeManagement.git

## Install from Source (Standard)

    git clone https://github.com/evolu-tion/GenomeManagement.git
    cd GenomeManagement
    pip install .

## Install for Development (Editable)

    pip install -e .

# Usage: Promoter of genes retrieving

Then run python script by used command line on windows or unix:

    promoter-retrieve \
    	--output <output_file.fa> \
    	--output_format <fasta/gff> \
    	--genome <genome.fa> \
    	--gff <genome.gff> \
    	--type <TLS/TSS> \
    	--upstream <bp> \
    	--downstream <bp> \
    	--all_gene <Y/N> \
    	--selected_gene_list <gene_list.txt, is optional if all_gene is N> \
    	--remove_n_gap <Y/N> \
    	--min_length <default is 100 bp>

# Usage: Get protein or gene sequences from genome

    proteins-retrieve \
    	--input <input_fasta_file> \
    	--list_of_interest <list_of_protein_id.txt> \
    	--output <output_file.fa>

# Usage: PlantPAN3

    plantpan3

# Usage: Get genome statistic

    get-genome-statistic \
    	--genome <genome.fa>

# License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
