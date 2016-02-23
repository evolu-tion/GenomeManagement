#!/usr/bin/env python
""" 
Copyright (c) 2016 King Mongkut's University technology Thonburi
Author: Nattawet Sriwichai
Contact: nattawet.sri@mail.kmutt.ac.th
License: MIT License

The MIT License

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
"""

############################### Initial Configuration ###########################################
# Genome input files including: 
# 1) Genome sequence file (fasta format)
# 2) Gene feature format  (gff3 format)
# 3) List of genes (tab-delimited format) is optional when you want to custome list
# cassava genome
file_genome_seq = 'Mesculenta_305_v6.fa'
file_gff = 'Mesculenta_305_v6.1.gene.gff3'
file_custom_list_of_gene = 'list_gene.txt'

# Configuration promoter properties
# 1) Promoter stating form 'TSS' or 'TLS' 
# 2) Upstream length form 'TSS' or 'TLS' 
# 3) Downstream length form 'TSS' or 'TLS' 
# 4) If you custom list of interested gene please type 'No', wherease if you want to retrieved whole genome please type 'Yes' 

start_promoter_from = 'TSS'
upstream = 2000
downstream = 0
promoter_minimum_length = 500
all_promoter_in_genome = 'Yes'

# Output file including:
# 1) Promoter sequences
# 2) List of genes cannot retrieved promoters
output_file_promoter = 'out/promoter.fa'
output_file_list_no_promoter = 'out/list_no_promoter.txt'

##################################################################################################



##################################### Python Main Program ########################################

import os
import sys
from seq_manage import Fasta_manager
from seq_manage import Gff_manager
from seq_manage import Genome_manager


def main():
	if not os.path.exists(file_genome_seq):
		print("Location of genome sequence file is not correct")
		exit()
	if not os.path.exists(file_gff):
		print("Location of gene feature file is not correct")
		exit()
	if all_promoter_in_genome == 'No' and not os.path.exists(file_custom_list_of_gene):
		print("Location of list of genes is not correct")
		exit()

	os.makedirs(os.path.dirname(output_file_promoter), exist_ok=True)
	os.makedirs(os.path.dirname(output_file_list_no_promoter), exist_ok=True)
	out_promoter = open(output_file_promoter, 'w')

	if(all_promoter_in_genome == 'No'):
		custom_gene_list = open(file_custom_list_of_gene).read().splitlines()
		
		# Check gene list in genomes
		print("Checking gene in", file_custom_list_of_gene)
		gff = Gff_manager(file_gff)
		for gene_name in custom_gene_list:
			if gff.checkGene(gene_name) == False:
				print(gene_name, "not in genome, please checking before new run.")
				exit()

		# Read Genome sequence and annotation
		print("Read Genome files..")
		genome = Genome_manager(file_genome_seq, file_gff)

		# FInding promoter on custom list
		if(start_promoter_from == 'TLS'):
			for gene_name in custom_gene_list:
				print("Finding promoter from TLS of ", gene_name)
				out_promoter.write(genome.getPromoterOfGeneFromTLS(gene_name, upstream, downstream, promoter_minimum_length))
		elif(start_promoter_from == 'TSS'):
			for gene_name in custom_gene_list:
				print("Finding promoter from TSS  of ", gene_name)
				out_promoter.write(genome.getPromoterOfGeneFromTSS(gene_name, upstream, downstream, promoter_minimum_length))
	else:

		# Finding promoter of all genes in genome
		print("Read Genome files..")
		genome = Genome_manager(file_genome_seq, file_gff)
		
		if(start_promoter_from == 'TLS'):
			for gene_name in genome.getGeneList():
				print("Finding promoter from TLS  of ", gene_name)
				out_promoter.write(genome.getPromoterOfGeneFromTLS(gene_name, upstream, downstream, promoter_minimum_length))
		elif(start_promoter_from == 'TSS'):
			for gene_name in genome.getGeneList():
				print("Finding promoter from TSS  of ", gene_name)
				out_promoter.write(genome.getPromoterOfGeneFromTSS(gene_name, upstream, downstream, promoter_minimum_length))

	print("Write list of gene no promoter")
	error_promoter = open(output_file_list_no_promoter, 'w')
	for gene_name in genome.getListOfGeneNoPromoter():
		error_promoter.write(gene_name +'\n')

if __name__ == "__main__":
    main()