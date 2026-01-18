#!/usr/bin/python




version = "GenomeManagement_v2.0.0"
import os
import sys
from optparse import OptionParser
from .seq_manage import FastaManager


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

def main():
    if not options.file_output:
        sys.exit("Missing output file, -o <FILE> or --output=<FILE>")
    if not options.file_protein_seq or not os.path.exists(options.file_protein_seq):
        sys.exit("Missing FASTA of reference protein sequence file, -i <FILE> or --input=<FILE>")
    if not options.list_of_interest:
        sys.exit("Missing retrieving protein list file, -l <FILE> or --list_of_interest=<FILE>")

    list_of_interest = options.list_of_interest
    file_output = options.file_output
    file_protein_seq = options.file_protein_seq

    protein = FastaManager(file_protein_seq)
    
    with open(list_of_interest, 'r') as custom_gene_list, open(file_output, 'w') as output_file:
        for each_gene in custom_gene_list:
            parts = each_gene.split()
            if len(parts) > 1:
                protein_id = parts[0]
                protein_symbol = parts[1]
            else:
                protein_id = parts[0]

            if protein.check_chromosome(protein_id):
                print("Found protein:", protein_id)
                sequence = protein.get_chr_sequence(protein_id)
                if len(parts) > 1:
                    output_file.write(f">{protein_id} {protein_symbol}\n{sequence}\n")
                else:
                    output_file.write(f">{protein_id}\n{sequence}\n")
            else:
                print(protein_id, "did not found protein reference file")

if __name__ == "__main__":
    main()
