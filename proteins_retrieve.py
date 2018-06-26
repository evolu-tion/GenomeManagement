#!/usr/bin/python

""" 
Copyright (c) 2016 King Mongkut's University technology Thonburi
Author: Nattawet Sriwichai
Contact: nattawet.sri@mail.kmutt.ac.th
Version: 1.3a 2017-05-31
License: MIT License
"""

version = "GenomeManagement_v1.3"
import os
import sys
from optparse import OptionParser
from seq_manage import Fasta_manager


if "--version" in sys.argv[1:] or "-v" in sys.argv[1:]:
	# TODO - Capture version of Select_representative_miRNA
	print(version)
	sys.exit(0)

# Parse Command Line
usage = """

Description:
This script designed for retrieving promoter sequences.

usage:
$ python proteins_retrieve.py \\
	--input <protein.fa> \\
	--list_of_interest <protein_list.txt> \\
	--output <output_file> \\
"""

parser = OptionParser(usage=usage)
parser.add_option("-o", "--output", dest="file_output",
	default=None, metavar="FILE",
	help="Output file name")
parser.add_option("-i", "--input", dest="file_protein_seq",
	default=None, metavar="FILE",
	help="Input FASTA of protein sequence file")
parser.add_option("-l", "--list_of_interest", dest="list_of_interest",
	default=None, metavar="FILE",
	help="List of selecting genes for retrieving promoter in text file (is optional if do not selecting all genome)")
options,args = parser.parse_args()

if not options.file_output:
	sys.exit("Missing output file, -o <FILE> or --output=<FILE>")
if not options.file_protein_seq or not os.path.exists(options.file_protein_seq):
	sys.exit("Missing FASTA of reference protein sequence file, -i <FILE> or --input=<FILE>")
if not options.list_of_interest:
	sys.exit("Missing retrieving protein list file, -l <FILE> or --list_of_interest=<FILE>")

list_of_interest = options.list_of_interest
file_output = options.file_output
file_protein_seq = options.file_protein_seq

# Output file including:
# 1) Protein sequence on list of protein

##################################### Python Main Program ########################################

def main():
	protein = Fasta_manager(file_protein_seq)
	custom_gene_list = open(list_of_interest, 'r')
	output_file = open(file_output, 'w')

	# FInding promoter on custom list
	for each_gene in custom_gene_list:
		print(each_gene)
		if len(each_gene.split()) > 1:
			protein_id = each_gene.split()[0]
			protein_symbol = each_gene.split()[1]
		else:
			protein_id = each_gene

		if protein.checkChromosome(protein_id):
			print("Found protein:", protein_id)
			sequence = protein.getChrSequence(protein_id)
			# print(">" + protein_id + "\n" + sequence + "\n")
			if len(each_gene.split()) > 1:
				output_file.write(">" + protein_id + " "+ protein_symbol + "\n" + sequence + "\n")
			else:
				output_file.write(">" + protein_id + "\n" + sequence + "\n")
		else:
			print(protein_id, "did not found protein reference file")
	output_file.close()

main()
