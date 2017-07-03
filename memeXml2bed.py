#!/usr/bin/python3

""" 
Copyright (c) 2016 King Mongkut's University technology Thonburi
Author: Nattawet Sriwichai
Contact: nattawet.sri@mail.kmutt.ac.th
Version: 1.0
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

import os
import sys
from optparse import OptionParser

if "--version" in sys.argv[1:]:
	# TODO - Capture version of Select_representative_miRNA
	print("memeXml2bed version 1.0")
	sys.exit(0)

# Parse Command Line
usage = """

Description:
This script designed for convert MEME output to BED.

usage:
$ python3 memeXml2bed.py \\
	--output <out.bed> \\
	--input <meme.xml> \\
	--pvalue <1e-2> \\
"""

parser = OptionParser(usage=usage)
parser.add_option("-o", "--output", dest="file_output",
	default='meme.bed', metavar="FILE",
	help="Output file name of BED, default is meme.out")
parser.add_option("-i", "--input", dest="file_xml",
	default='meme.xml', metavar="FILE",
	help="Input XML file from MEME, default is meme.xml")
parser.add_option("-p", "--pvalue", dest="pvalue",
	default=1e-2, metavar="FILE",
	help="Cut off motif scan p-value that less than, default is 1e-2")

options,args = parser.parse_args()
file_output = options.file_output
file_xml = options.file_xml
pvalue = float(options.pvalue)

if not options.file_output:
	sys.exit("Missing output BED file, -o <FILE> or --output=<FILE>")
if not options.file_xml or not os.path.exists(options.file_xml):
	sys.exit("Missing input XML file from MEME, -i <FILE> or --input=<FILE>")


import xml.etree.ElementTree as ET
tree = ET.parse(file_xml)
root = tree.getroot()

# append input sequence into dict
sequence = {}
for i in range(1, len(root[0])-1):
	sequence[root[0][i].attrib['id']] = root[0][i].attrib

# read command line
command_line = root[1][0].text

# append motif information into dict
motif = {}
for i in range(0, len(root[2])):
	motif[root[2][i].attrib['id']] = root[2][i].attrib

# append motif_sites into dict
motif_sites = {}
for i in range(0, len(root[3])):
	motif_sites_in_sequence = []
	for j in range(len(root[3][i])):
		motif_sites_in_sequence.append(root[3][i][j].attrib)
	sequence[root[3][i].attrib['sequence_id']]['motif_sites'] = motif_sites_in_sequence

# print motif location

file_out = open(file_output, 'w')

for i in sequence:
	name = sequence[i]['name']
	name = name.split('|')[0]
	sequence_length = int(sequence[i]['length'])

	file_out.write(name + '\t' + str(sequence_length) + '\t' + str(sequence_length+10)  + '\t' + 'CDS\n')
	# print(name, sequence_length, sequence_length+10, 'CDS', sep='\t')
	for j in range(len(sequence[i]['motif_sites'])):
		motif_id = sequence[i]['motif_sites'][j]['motif_id']
		motif_position = sequence[i]['motif_sites'][j]['position']
		motif_pvalue = float(sequence[i]['motif_sites'][j]['pvalue'])
		motif_length = motif[motif_id]['width']
		if motif_pvalue < pvalue:
			file_out.write(name + '\t' + str(int(motif_position) + 1) + '\t' + str(int(motif_position) + int(motif_length)) + '\t' + motif_id + '\n')
			# print(name, int(motif_position) + 1, int(motif_position) + int(motif_length), motif_id, sep='\t')
