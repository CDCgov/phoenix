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
# v2.0 (4/10/2023)
#
# Originally, created by Erisa Sula (nvd4@cdc.gov) with v2 written by Jill Hagey (qpk9@cdc.gov)
#

from Bio import SeqIO
import sys
import os
import argparse
import re

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to rename contigs in assemblies')
	parser.add_argument('-i', '--input', required=True, help='input fasta filename')
	parser.add_argument('-o', '--output', required=True, help='output filename')
	parser.add_argument('-n', '--name', dest="name", required=True, help='filename')
	parser.add_argument('--reverse', help='returns formatted header to original', action='store_true')
	return parser.parse_args()

args=parseArgs()

def spades_rename(input_file, name, reverse, output):
	"""
	Typical spades output: >NODE_1_length_600507_cov_51.436379
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#_depth_depthx)
	"""
	sequences = []
	if not reverse:
		for record in SeqIO.parse(input_file,"fasta"):
			record.id = record.id.split("_cov")[0].replace("NODE",name)
			record.description = ""
			sequences.append(record)
		SeqIO.write(sequences, output, "fasta")
	else:
		for record in SeqIO.parse(input_file,"fasta"):
			record.id = record.id.replace(name,"NODE")+"_cov_X"
			record.description = ""
			sequences.append(record)
		SeqIO.write(sequences, output, "fasta")

def shovill_rename(input_file, name, output):
	"""
	Typical shovill output: >contig00001 len=1423313 cov=21.0 corr=0 origname=NODE_1_length_1423313_cov_21.008026_pilon sw=shovill-spades/1.1.0 date=20230327
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#_depth_depthx)
	"""
	sequences = []
	for record in SeqIO.parse(input_file,"fasta"):
		record.id = re.search('origname=(.*) sw=', record.description).group(1)
		record.id = record.id.split("_cov")[0].replace("NODE",name)
		record.description = ""
		sequences.append(record)
	SeqIO.write(sequences, output, "fasta")

def skesa_rename(input_file, name, output):
	"""
	Typical skesa output: >Contig_1_42.0537_Circ [topology=circular]
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#_depth_depthx)
	"""
	sequences = []
	for record in SeqIO.parse(input_file,"fasta"):
		record.id = record.id.split("_")[0] + "_" + record.id.split("_")[1] + "_length_" + str(len(record.seq))
		record.id = record.id.replace("Contig",name)
		record.description = ""
		sequences.append(record)
	SeqIO.write(sequences, output, "fasta")

def detect_assemblier(input_file, name):
	for record in SeqIO.parse(input_file,"fasta"):
		if record.id.startswith("contig00001"):
			assemblier = "shovill"
		elif record.id.startswith("NODE_1"):
			assemblier = "spades"
		elif record.id.startswith("Contig_1"):
			assemblier = "skesa"
		elif record.id.startswith(name + "_1_length_"):
			assemblier = "correct name"
	return assemblier

def rename_file(input_file, output):
	sequences = []
	for record in SeqIO.parse(input_file,"fasta"):
		sequences.append(record)
	SeqIO.write(sequences, output, "fasta")

def main():
	args = parseArgs()
	assemblier = detect_assemblier(args.input, args.name)
	if assemblier == "spades":
		spades_rename(args.input, args.name, args.reverse, args.output)
	elif assemblier == "skesa":
		skesa_rename(args.input, args.name, args.output)
	elif assemblier == "shovill":
		shovill_rename(args.input, args.name, args.output)
	elif assemblier == "correct name": # if the name looks like the correct scheme then just rename and move on
		rename_file(args.input, args.output)
	else:
		print("ERROR: The assemblier used could not be determined. PHoeNIx supports assemblies from either Skesa, Shovill or SPAdes.")

if __name__ == '__main__':
	main()