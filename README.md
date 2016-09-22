# Genome Management Package
This tool is copyright 2016 by Nattawet Sriwichai, King Mongkut's University Technology Thonburi (Bioinformatics and Systems Biology), Th. All rights reserved. See the licence text below.

# Usage: Promoter of gene retrieving
Edit the configuration on ``promoter_retrieve.py`` of program required and set parameter including

1. Genome File: Genome sequence file, Gene feature file, List of genes to retreive promoter

2. Promoter Properties: the stating point of promoter is 'TSS' or 'TLS', length of upstream and downstream.

Then run python script by used command line on windows or unix::

    python promoter_retrieve.py \\
		--out_promoter <output_file.fa> \\
		--genome <genome.fa> \\
		--gff <genome.gff> \\
		--type <TLS/TSS> \\
		--upstream <bp> \\
		--downstream <bp> \\
		--all_gene <Y/N> \\
		--list_of_interest <gene_list.txt> \\
		--remove_n_gap <Y/N>

# Licence (MIT)
Copyright (c) 2016 King Mongkut's University technology Thonburi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE. 
