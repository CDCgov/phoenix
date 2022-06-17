#!/usr/bin/env python3

#
# Description: Changes headers in SPAdes assembly fasta from contig# length=length# depth=depthx to Name_contig#_length_length#_depth_depthx
#
# Usage: ./rename_fasta_headers.py -i input.fasta -o output.fasta
#
# Output location: parameter
#
# Modules required: Biopython must be available in python instance
#
# v1.0 (10/3/2019)
#
# Created by Erisa Sula (nvd4@cdc.gov)
#

from Bio import SeqIO
import sys
import os
import argparse

#print("Starting")
#Create an arg parser...someday
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to rename contigs in assemblies')
	parser.add_argument('-i', '--input', required=True, help='input fasta filename')
	parser.add_argument('-o', '--output', required=True, help='output filename')
	parser.add_argument('-n', '--name', dest="name", required=True, help='filename')
	parser.add_argument('--reverse', help='returns formatted header to original', action='store_true')
	return parser.parse_args()

args=parseArgs()
sequences = []

if not args.reverse:
	for record in SeqIO.parse(args.input,"fasta"):
		#print(record.id)
		#print(name)
		record.id = record.id.split("_cov")[0].replace("NODE",args.name)
		#print(record.id)
		record.description = ""
	#    print(record.description)
	#    print(record)
		sequences.append(record)

	SeqIO.write(sequences, args.output, "fasta")
else:
	#print("REVERSE")
	#print(name)
	for record in SeqIO.parse(args.input,"fasta"):
		#print(record.id)
		#print(name)
		record.id = record.id.replace(args.name,"NODE")+"_cov_X"
		#print(record.id)
		record.description = ""
		#print(record.description)
		#print(record)
		sequences.append(record)

	SeqIO.write(sequences, args.output, "fasta")
