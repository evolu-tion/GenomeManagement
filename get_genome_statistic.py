#!/usr/bin/env python

import os
import sys
from optparse import OptionParser

if "--version" in sys.argv[1:]:
	# TODO - Capture version of get genome statistic
	print("get_genome_statistic v1")
	sys.exit(0)

# Parse Command Line
usage = """

Description:
This script designed for geting genome statistic.

usage:
$ python get_genome_statistic.py \\
	--genome <genome.fa> \\
"""

parser = OptionParser(usage=usage)
parser.add_option("-g", "--genome", dest="file_genome_seq",
	default=None, metavar="FILE",
	help="Input FASTA of genome sequence file")

options,args = parser.parse_args()

if not options.file_genome_seq or not os.path.exists(options.file_genome_seq):
	sys.exit("Missing FASTA of genome sequence file, -g <FILE> or --genome=<FILE>")


############################### Initial Configuration ###########################################
# Genome input files including: 
# 1) Genome sequence file (fasta format)
file_genome_seq = options.file_genome_seq

##################################################################################################



##################################### Python Main Program ########################################

import os
import sys
from seq_manage import Fasta_manager


def main():
	genome = Fasta_manager(file_genome_seq, True)

if __name__ == "__main__":
    main()
