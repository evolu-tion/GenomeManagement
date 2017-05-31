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


file_custom_list_of_protein = "genome/list_of_proteins.txt"

# Output file including:
# 1) Protein sequence on list of protein
output_file_protein_file = "Athaliana_AGPase.fa"


##################################### Python Main Program ########################################
import os
import sys
from seq_manage import Fasta_manager

def main():
	protein = Fasta_manager("genome/Athaliana_167_TAIR10.protein.fa")
	custom_gene_list = open(file_custom_list_of_protein).read().splitlines()
	output_file = open("output_file_protein_file", 'w')

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