#!/usr/bin/env python

import os
import sys
from optparse import OptionParser
from .seq_manage import FastaManager

if "--version" in sys.argv[1:]:
	# TODO - Capture version of get genome statistic
	print("get_genome_statistic v2.0.0")
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


def main():
    if not options.file_genome_seq or not os.path.exists(options.file_genome_seq):
        sys.exit("Missing FASTA of genome sequence file, -g <FILE> or --genome=<FILE>")
    
    file_genome_seq = options.file_genome_seq
    genome = FastaManager(file_genome_seq, True)

if __name__ == "__main__":
    main()
