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
import os
import sys
from seq_manage import Fasta_manager

def main():
	protein = Fasta_manager("genome/Athaliana_167_TAIR10.protein.fa")
	custom_gene_list = open(list_of_interest).read().splitlines()
	output_file = open(file_output, 'w')

	# FInding promoter on custom list
	for protein_id in custom_gene_list:
		if protein.checkChromosome(protein_id):
			print("Found protein:", protein_id)
			sequence = protein.getSequence(protein_id)
			# print(">" + protein_id + "\n" + sequence + "\n")
			output_file.write(">" + protein_id + "\n" + sequence + "\n")
		else:
			print(protein_id, "did not found protein reference file")

if __name__ == "__main__":
    main()