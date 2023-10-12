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
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#)
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
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#)
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
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#)
	"""
	sequences = []
	for record in SeqIO.parse(input_file,"fasta"):
		record.id = record.id.split("_")[0] + "_" + record.id.split("_")[1] + "_length_" + str(len(record.seq))
		record.id = record.id.replace("Contig",name)
		record.description = ""
		sequences.append(record)
	SeqIO.write(sequences, output, "fasta")

def unicycler_rename(input_file, name, output):
	"""
	Typical unicycler output: >1 length=2836348 depth=1.00x circular=tru
	Convert to: >1921706_1_length_2836348 (Name_contig#_length_length#_depth_depthx)
	"""
	sequences = []
	for record in SeqIO.parse(input_file,"fasta"):
		record.id = name + "_" + str(record.id) + "_length_" + str(record.description.split(" ")[1].replace("length=", ""))
		record.description = ""
		sequences.append(record)
	SeqIO.write(sequences, output, "fasta")

def flye_rename(input_file, name, output):
	"""
	Typical Flye output: >contig_1
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#)
	"""
	sequences = []
	for record in SeqIO.parse(input_file,"fasta"):
		record.id = name + "_" + record.id.split("_")[1] + "_length_" + str(len(record.seq))
		#record.id = record.id.replace("contig",name)
		record.description = ""
		sequences.append(record)
	SeqIO.write(sequences, output, "fasta")

def trycycler_rename(input_file, name, output):
	"""
	Typical Flye output: >cluster_001_consensus_polypolish or >cluster_001_consensus
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#)
	"""
	sequences = []
	for record in SeqIO.parse(input_file,"fasta"):
		record.id = name + "_" + str(int(record.id.split("_")[1])) + "_length_" + str(len(record.seq))
		record.id = record.id.replace("Contig",name)
		record.description = ""
		sequences.append(record)
	SeqIO.write(sequences, output, "fasta")

def unknown_rename(input_file, output, name):
	"""
	Convert to: >1921706_1_length_600507 (Name_contig#_length_length#)
	"""
	sequences = []
	count = 0
	for record in SeqIO.parse(input_file,"fasta"):
		record.id = name + "_" + str(count) + "_length_" + str(len(record.seq))
		record.description = ""
		sequences.append(record)
		count = count + 1
	SeqIO.write(sequences, output, "fasta")

def detect_assemblier(input_file, name):
	for record in SeqIO.parse(input_file,"fasta"):
		if record.id.startswith("contig00001"):
			assembler = "shovill"
		elif record.id.startswith("NODE_1"):
			assembler = "spades"
		elif record.id.startswith("Contig_1"):
			assembler = "skesa"
		elif record.id.startswith("contig_1"):
			assembler = "flye"
		elif record.id.startswith("cluster_001"):
			assembler = "trycycler"
		elif record.id.startswith("1") and record.description.startswith("1 length="):
			assembler = "unicycler"
		elif record.id.startswith(name + "_1_length_"):
			assembler = "correct name"
		else: # if none of these then leave blank for error to be reported
			assembler = "Unknown"
		break # we only need to do this for the first record
	return assembler

def rename_file(input_file, output):
	sequences = []
	for record in SeqIO.parse(input_file,"fasta"):
		sequences.append(record)
	SeqIO.write(sequences, output, "fasta")

def main():
	args = parseArgs()
	assembler = detect_assemblier(args.input, args.name)
	if assembler == "spades":
		spades_rename(args.input, args.name, args.reverse, args.output)
	elif assembler == "skesa":
		skesa_rename(args.input, args.name, args.output)
	elif assembler == "shovill":
		shovill_rename(args.input, args.name, args.output)
	elif assembler == "unicycler":
		unicycler_rename(args.input, args.name, args.output)
	elif assembler == "flye": 
		flye_rename(args.input, args.name, args.output)
	elif assembler == "trycycler":
		trycycler_rename(args.input, args.name, args.output)
	elif assembler == "correct name": # if the name looks like the correct scheme then just rename and move on
		rename_file(args.input, args.output)
	elif assembler == "Unknown": # if the name isn't recogized then
		unknown_rename(args.input, args.output, args.name)
	else:
		print("ERROR: The assembler used could not be determined. PHoeNIx supports assemblies from either Skesa, Shovill, SPAdes, Unicycler,Trycycler and Flye. If you used one of these open a github issue to report the problem.")

if __name__ == '__main__':
	main()