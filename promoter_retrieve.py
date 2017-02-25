#!/usr/bin/python

""" 
Copyright (c) 2016 King Mongkut's University technology Thonburi
Author: Nattawet Sriwichai
Contact: nattawet.sri@mail.kmutt.ac.th
Version: 1.3 2017-02-25
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
version = "GenomeManagement_v1.3"

import os
import sys
from optparse import OptionParser

if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
	# TODO - Capture version of Select_representative_miRNA
	print(version)
	sys.exit(0)

# Parse Command Line
usage = """

Description:
This script designed for retrieving promoter sequences.

usage:
$ python promoter_retrieve.py \\
	--output <output_file> \\
	--output_format <fasta/gff> \\
	--genome <genome.fa> \\
	--gff <genome.gff> \\
	--type <TLS/TSS> \\
	--upstream <bp> \\
	--downstream <bp> \\
	--all_gene <Y/N> \\
	--list_of_interest <gene_list.txt> \\
	--remove_n_gap <Y/N>
"""

parser = OptionParser(usage=usage)
parser.add_option("-o", "--output", dest="file_output",
	default=None, metavar="FILE",
	help="Output file name")
parser.add_option("-p", "--output_format", dest="output_format",
	default=None, metavar="<fasta/fa/FASTA/GFF/gff/GFF3/gff3>",
	help="Output file format")
parser.add_option("-g", "--genome", dest="file_genome_seq",
	default=None, metavar="FILE",
	help="Input FASTA of genome sequence file")
parser.add_option("-f", "--gff", dest="file_gff",
	default=None, metavar="FILE",
	help="Input GFF annotation file")
parser.add_option("-t", "--type", dest="type",
	default=None, metavar="<TSS/TLS>",
	help="Selecting the promoter sequence from TSS or TLS")
parser.add_option("-u", "--upstream", dest="upstream",
	default=1500, metavar="int",
	help="Length of upstream from TLS or TSS (bps)")
parser.add_option("-d", "--downstream", dest="downstream",
	default=0, metavar="int",
	help="Length of downstream from TLS or TSS (bps)")
parser.add_option("-a", "--all_gene", dest="all_gene",
	default="Yes", metavar="<Y/N>",
	help="If you want to retrieve all genes in genome please type <Yes> or <No>")
parser.add_option("-l", "--selected_gene_list", dest="gene_list",
	default=None, metavar="FILE",
	help="List of selecting genes for retrieving promoter in text file (is optional if do not selecting all genome)")
parser.add_option("-r", "--remove_n_gap", dest="remove_N_gap",
	default='N', metavar="<Y/N>",
	help="Remove N gap (Y/N)")
parser.add_option("-m", "--min_length", dest="min_length",
	default=100, metavar="int",
	help="Minimum number of promoter length (bp), default=100")
options,args = parser.parse_args()

if not options.file_output:
	sys.exit("Missing output file, -o <FILE> or --output=<FILE>")
if not options.output_format:
	sys.exit("Missing output format file, -p <fasta/gff> or --output_format=<fasta/gff>")
if not options.file_genome_seq or not os.path.exists(options.file_genome_seq):
	sys.exit("Missing FASTA of genome sequence file, -g <FILE> or --genome=<FILE>")
if not options.file_gff or not os.path.exists(options.file_gff):
	sys.exit("Missing GFF3 annotation file of genome, -f <FILE> or --gff=<FILE>")
if not options.type:
	sys.exit("Missing information of promoter sequence from TSS or TLS, -t <TLS/TSS> or --type=<TLS/TSS>")
if not options.upstream:
	sys.exit("Missing length of upstream from TLS or TSS (bps), -u int or --upstream=int")
if not options.downstream:
	sys.exit("Missing length of downstream from TLS or TSS (bps), -d int or --downstream=int")
if not options.all_gene:
	sys.exit("Missing retrieve all genes in genome or not, -a int or --all_gene=<Yes/No>")

############################### Initial Configuration ###########################################
# Genome input files including: 
# 1) Genome sequence file (fasta format)
# 2) Gene feature format  (gff3 format)
# 3) List of genes (tab-delimited format) is optional when you want to custome list
# cassava genome
file_genome_seq = options.file_genome_seq
file_gff = options.file_gff
file_custom_list_of_gene = options.gene_list

# Configuration promoter properties
# 1) Promoter stating form 'TSS' or 'TLS' 
# 2) Upstream length form 'TSS' or 'TLS' 
# 3) Downstream length form 'TSS' or 'TLS' 
# 4) If you custom list of interested gene please type 'No', wherease if you want to retrieved whole genome please type 'Yes' 
# 5) If you want to remove N gap in begin sequence please type 'Yes', wherease type 'No'

start_promoter_from = options.type
upstream = int(options.upstream)
downstream = int(options.downstream)
promoter_minimum_length = 100
all_promoter_in_genome = options.gene_list
removed_N_gap_begins_promoter = options.remove_N_gap

# Output file including:
# 1) Promoter sequences
# 2) List of genes cannot retrieved promoters
output_file_promoter = options.file_output
output_format = options.output_format
output_file_list_no_promoter = options.file_output + '_list_no_promoter.txt'

##################################################################################################



##################################### Python Main Program ########################################

import os
import sys
from seq_manage import Fasta_manager
from seq_manage import Gff_manager
from seq_manage import Genome_manager


def main():
	if (options.all_gene.upper() == 'N' or options.all_gene.upper == 'NO') and not os.path.exists(file_custom_list_of_gene):
		print("Location of list of genes is not correct")
		exit()
	
	if removed_N_gap_begins_promoter == 'Y' or removed_N_gap_begins_promoter == 'y':
		removed_N_gap = True
	elif removed_N_gap_begins_promoter == 'N' or removed_N_gap_begins_promoter == 'n':
		removed_N_gap = False
	else:
		print("Missing --remove_N_gap=<Y/N>")
		exit()


	os.makedirs(os.path.dirname(output_file_promoter), exist_ok=True)
	os.makedirs(os.path.dirname(output_file_list_no_promoter), exist_ok=True)
	out_promoter = open(output_file_promoter, 'w')

	if(options.all_gene.upper() == 'N' or options.all_gene.upper == 'NO'):
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
				out_promoter.write(genome.getPromoterOfGeneFromTLS(gene_name, upstream, downstream, promoter_minimum_length, removed_N_gap, output_format))
		elif(start_promoter_from == 'TSS'):
			for gene_name in custom_gene_list:
				print("Finding promoter from TSS  of ", gene_name)
				out_promoter.write(genome.getPromoterOfGeneFromTSS(gene_name, upstream, downstream, promoter_minimum_length, removed_N_gap, output_format))
	else:

		# Finding promoter of all genes in genome
		print("Read Genome files..")
		genome = Genome_manager(file_genome_seq, file_gff)
		
		if(start_promoter_from == 'TLS'):
			for gene_name in genome.getGeneList():
				print("Finding promoter from TLS  of ", gene_name)
				out_promoter.write(genome.getPromoterOfGeneFromTLS(gene_name, upstream, downstream, promoter_minimum_length, removed_N_gap, output_format))
		elif(start_promoter_from == 'TSS'):
			for gene_name in genome.getGeneList():
				print("Finding promoter from TSS  of ", gene_name)
				out_promoter.write(genome.getPromoterOfGeneFromTSS(gene_name, upstream, downstream, promoter_minimum_length, removed_N_gap, output_format))

	print("Write list of gene no promoter")
	error_promoter = open(output_file_list_no_promoter, 'w')
	for gene_name in genome.getListOfGeneNoPromoter():
		error_promoter.write(gene_name +'\n')

if __name__ == "__main__":
    main()
