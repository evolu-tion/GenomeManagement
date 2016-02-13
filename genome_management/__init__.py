# Author: Nattawet Sriwichai
# Author-email: nattawet.sri@mail.kmutt.ac.th
# https://github.com/evolu-tion/GenomeManagement
# Copyright (c) 2016 Nattawet Sriwichai
# License: LICENSE.txt 

############################### Initial Configuration ###########################################
# Genome input files including: 
# 1) Genome sequence file (fasta format)
# 2) Gene feature format  (gff3 format)
# 3) List of genes (tab-delimited format) is optional when you want to custome list
# 4) Output file of promoter (fasta format)
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
all_promoter_in_genome = 'All'

# Output file including:
# 1) Promoter sequences
# 2) List of genes cannot retrieved promoters
output_file_promoter = 'out/promoter.fa'
output_file_list_no_promoter = 'out/list_no_promoter.txt'

##################################################################################################



##################################### Python Main Program ########################################


from seq_manage import Fasta_manager
from seq_manage import Gff_manager
from seq_manage import Genome_manager

def main():
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