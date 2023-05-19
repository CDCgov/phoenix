#!/usr/bin/env python3

#
# Description: Script to go through an isolates MLST output and check and reformat any instances of partial or multiple profiles found. SUB for profile or allele means that there is
#	   a novel allele or the profile has not been assigned yet, hence something needs to be submitted to pubmlst. AU (Allele Unknown) implies that an allele can not be determined
#	   no profile can be determined from the current assembly/qc_reads
#
# Usage: is python3 ./fix_MLST.py [-i input_torsten_MLST_file] [-s input_srst2_mlst_file] -t taxonomy_file
#
# Modules required: None
#
# v3 (05/08/2023)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#from asyncio.windows_events import NULL
from locale import currency
import sys
#disable cache usage in the Python so __pycache__ isn't formed. If you don't do this using 'nextflow run cdcgov/phoenix...' a second time will causes and error
sys.dont_write_bytecode = True
import os
#import glob
#import math
import copy
import argparse
import collections
#import itertools as it
#from pathlib import Path
#import local_MLST_converter
import xml.dom.minidom as xml
#import urllib as url
#import re
#from urlparse import urlparse # this is for python2 line below is updated package for python3, Python3 needed for terra
from urllib.parse import urlparse
#import urllib.request
from datetime import datetime
from os import path
#import subprocess
#import ssl

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to check MLST types for duplicate alleles and implications on final typing')
	parser.add_argument('-i', '--input', required=False, help='input mlst filename')
	parser.add_argument('-s', '--srst2', required=False, help='srst2 input file')
	parser.add_argument('-t', '--taxonomy', required=True, help='Location of taxonomy file to pull right scheme')
	parser.add_argument('-d', '--mlst_database', required=True, help='Path to mlst db of db/pubmlst/schemes/ format')
	return parser.parse_args()


# main function that looks if all MLST types are defined for an outptu mlst file
def do_MLST_check(input_MLST_line_tuples, taxonomy_file, mlst_db_path):
	#location="/".join(taxonomy_file.split("/")[0:-1])+"/mlst"
	#print taxonomy_file
	isolate_name = taxonomy_file.split(".")[:-1]
	isolate_name = ".".join(isolate_name)
	tax_file = open(taxonomy_file, 'r')
	# Changing from real-time to set pull date
	#today=datetime.today().strftime('%Y-%m-%d')
	with open(mlst_db_path+'/db_version') as f:
		pull_date = f.readline().strip('\n')
	for line in tax_file:
		units=line.replace("\n", "").replace("\t","|").split("|")
		#print(len(units), units)
		if len(units) == 1:
			continue
		#print(line)
		tier=units[0]
		info=units[1]
		if tier == "G:":
			genus = info
		elif tier == "s:":
			species = info
			print("Taxonomy:",genus,species)
			break
	tax_file.close()

	print(input_MLST_line_tuples)

	original_schemes = []
	checked_schemes = []
	sample="None"
	db_name="None"
	original_type="Introvert"

	### Creates a list of all schemes found in all files from the MLST_file_list
	for input_MLST_line in input_MLST_line_tuples:
		#MLST_file = open(input_MLST_line[2],'r')
		MLST_line = input_MLST_line[0]
		MLST_filetype = input_MLST_line[1]
		original_line=MLST_line.strip()
		if original_line is not None and original_line != '':
			allele_list=[]
			expanded_allele_list=[]
			allele_names=[]
			print(original_line)
			original_items=original_line.split("\t")
			print("Array of original itmes","\n".join(original_items))

			if MLST_filetype == "mlst":
				sample=original_items[0]
				db_name=original_items[1]
				original_type=original_items[2]

				# Default list size in case it is Empty
				list_size=0

				# Create a list of pertinent info and make a list of 'sets' of alleles
				if len(original_items) > 4:
					for allele in range(3, len(original_items)):
	#					   print(":", original_items[allele])
						allele_Identifier=original_items[allele].split("(")[0]
						alleles=original_items[allele].split("(")[1].split(")")[0].split(",")
						if len(alleles) > 1:
							original_type="-"
						allele_names.append(allele_Identifier)
						allele_list.append(alleles)

					# Create a template to build Scheme array
					expanded_allele_list.append([sample,db_name, original_type, len(allele_names), allele_names, [], "standard", pull_date])

					# Test first parse
	#				   print("e:",expanded_allele_list)
	#				   print("o:", allele_list)

					# Expand allele/scheme list if any multiples are encountered to make one list entry per possible scheme
					for allele in allele_list:
						if len(allele) == 1:
							for i in expanded_allele_list:
								print(i, allele, allele[0])
								i[5].append(allele[0])
						else:
							temp_allele_list=copy.deepcopy(expanded_allele_list)
							print("More than one found for ", allele, "template=", temp_allele_list)
							for i in range(len(allele)):
								if i == 0:
	#								   print("0 - Going to add:", allele[i])
									for j in expanded_allele_list:
	#									   print("to - ", j)
										j[5].append(allele[i])
	#									   print("After adding - ", j)
								else:
									temp2_allele_list=copy.deepcopy(temp_allele_list)
	#								   print(temp2_allele_list)
	#								   print(i,"- Going to add:", allele[i])
									for j in temp2_allele_list:
	#									   print("to - ", j)
										j[5].append(allele[i])
	#									   print("After adding - ", j)
									expanded_allele_list.append(copy.deepcopy(temp2_allele_list)[0])
	#					   print(expanded_allele_list)
				else:
					#expanded_allele_list.append([sample,db_name, "-", 7, ["-","-","-","-","-","-","-"], ["-","-","-","-","-","-","-"], "standard"])
					expanded_allele_list.append([sample,db_name, "-", 1, ["-"], ["-"], "standard", pull_date])


			#allele_list=[['1'], ['3'], ['189','3'], ['2'], ['2'], ['96','107'], ['3']]
			elif MLST_filetype == "srst2":
				print(original_items,len(original_items))
				sample=original_items[0]
				db_name=original_items[1]
				if db_name == "No match found":
					continue
				original_type=original_items[2]
				if len(original_items) > 7:
					print("Has # items:", len(original_items))
					for allele in range(7, len(original_items)):
						#print(allele, original_items[allele])
						allele_Identifier=original_items[allele].split("(")[0]
						alleles=original_items[allele].split("(")[1].split(")")[0].split(",")
						if len(alleles) > 1:
							original_type="-"
						allele_names.append(allele_Identifier)
						allele_list.append(alleles)

					expanded_allele_list.append([sample,db_name, original_type, len(allele_names), allele_names, [], "srst2", pull_date])

					for allele in allele_list:
						if len(allele) == 1:
							for i in expanded_allele_list:
								print(i, allele, allele[0])
								i[5].append(allele[0])
						else:
							temp_allele_list=copy.deepcopy(expanded_allele_list)
							print("More than one found for ", allele, "template=", temp_allele_list)
							for i in range(len(allele)):
								if i == 0:
	#								   print("0 - Going to add:", allele[i])
									for j in expanded_allele_list:
	#									   print("to - ", j)
										j[5].append(allele[i])
	#									   print("After adding - ", j)
								else:
									temp2_allele_list=copy.deepcopy(temp_allele_list)
	#								   print(temp2_allele_list)
	#								   print(i,"- Going to add:", allele[i])
									for j in temp2_allele_list:
	#									   print("to - ", j)
										j[5].append(allele[i])
	#									   print("After adding - ", j)
									expanded_allele_list.append(copy.deepcopy(temp2_allele_list)[0])
						print(expanded_allele_list)
				else:
					expanded_allele_list.append([sample,db_name, "-", 1, ["-"], ["-"], "srst2", pull_date])

			else:
				print("Unknown MLST filetype, can not continue")
				exit()
			#original_line=MLST_file.readline().strip("	 ")
			original_schemes.append(expanded_allele_list)

		#MLST_file.close()

	### Shows the list of what schemes were found
	print("Schemes found:", len(original_schemes))
	catted_scheme_list=[]
	for oscheme in original_schemes:
		print(oscheme)
		for profile in oscheme:
			#print profile
			catted_scheme_list.append(profile)

	print("# of catted schemes found:", len(catted_scheme_list))
	for i in range(0, len(catted_scheme_list)):
		print(i, catted_scheme_list[i])

	dupes = []
	for i in range(0,len(catted_scheme_list)):
		for j in range(0,len(catted_scheme_list)):
			if i == j:
				continue
	#		   print("Checking:", catted_scheme_list[i][5], catted_scheme_list[j][5])
			#print catted_scheme_list[i],catted_scheme_list[j]
			if collections.Counter(catted_scheme_list[i][5]) == collections.Counter(catted_scheme_list[j][5]):
				print("SAME!!!!", i, j)
				if catted_scheme_list[min(i,j)][6] != catted_scheme_list[max(i,j)][6]:
					print("a")
					if catted_scheme_list[min(i,j)][6] == "standard":
						print("b")
						if catted_scheme_list[max(i,j)][6] == "srst2":
							print("c")
							#catted_scheme_list[min(i,j)][4] = catted_scheme_list[max(i,j)][4]
							catted_scheme_list[min(i,j)][6] = "standard/srst2"
							dupes.append(max(i,j))
						elif catted_scheme_list[max(i,j)][6] == "standard/srst2":
							print("d")
							dupes.append(min(i,j))
					elif catted_scheme_list[min(i,j)][6] == "srst2":
						print("e")
						if catted_scheme_list[max(i,j)][6] == "standard":
							print("f")
							#catted_scheme_list[max(i,j)][4] = catted_scheme_list[min(i,j)][4]
							catted_scheme_list[max(i,j)][6] = "standard/srst2"
							dupes.append(min(i,j))
						elif catted_scheme_list[max(i,j)][6] == "standard/srst2":
							print("g")
							dupes.append(min(i,j))
					elif catted_scheme_list[min(i,j)][6] == "standard/srst2":
						print("h")
						dupes.append(max(i,j))
					else:
						print("i")
						print("Should never have something that is not standard, srst2, or standard/srst2")
				else:
					if max(i,j) not in dupes:
						print("j")
						dupes.append(max(i,j))
					else:
						print("Index already Found")
				print(dupes)
			else:
	#			   print(catted_scheme_list[i][5], "does not equal", catted_scheme_list[i][5])
				continue
		dedupped_dupes=list(set(dupes))
	dedupped_dupes.sort(reverse=True)

	for k in dedupped_dupes:
		catted_scheme_list.pop(k)

	print("Trimmed catted: ", catted_scheme_list)

	### Begin checking list for completeness/errors
	if len(original_schemes) == 0:
		print("No profiles found")

	#######
		### Print out empty file
	#######

	else:
		for original_scheme in catted_scheme_list:
			current_type = original_scheme[2]
			current_alleles = original_scheme[5]
			print(current_type, current_alleles)
			check = False
			bad_types = ["*", "?", "NF", "~"]
			if any(ext in current_type for ext in bad_types):
				original_scheme[2] = "Novel_profile"
				print("Marking", original_scheme, "as needing investigation.")
				check = True
			elif current_type == "-":
				if original_scheme[3] > 1:
					original_scheme[2] = "Novel_profile"
					check = True
				else:
					original_scheme[2] = "-"
					check = False
			elif current_type == "failed":
				if original_scheme[3] > 1:
					original_scheme[2] = "Failed"
				check=False
			else:
				current_type = int(current_type)
				#print("Current MLST type:", current_type)
				#list_size=original_scheme[3]
				#print("Allele_names:", original_scheme[4])
				#print("Alleles_found:", original_scheme[5])
				check = False
				print("Adding ", original_scheme, "to check_schemes")
				print("B",checked_schemes)
				checked_schemes.append(original_scheme)
				print("A",checked_schemes)
			if check:
				print("About to look into", original_scheme)
				bad_alleles = 0
				lookup_allele_profile = True
				for allele in original_scheme[5]:
					# Definitions
					#	   MLST
					#	   ~ : full length novel allele
					#	   ? : partial match (>min_cov & > min_ID). Default min_cov = 10, Default min_ID=95%
					#	   - : Allele is missing
					#
					#	   srst2
					#	   * : Full length match with 1+ SNP. Novel
					#	   ? : edge depth is below N or average depth is below X. Default edge_depth = 2, Default average_depth = 5
					#	   - : No allele assigned, usually because no alleles achieved >90% coverage



					if '*' in allele or '?' in allele or '~' in allele or '-' in allele:
						original_scheme[2] = "Novel_allele"
						bad_alleles += 1
						lookup_allele_profile = False
				## Reapply dash showing no db or proximal scheme was ever found
				if collections.Counter(original_scheme[5]) == collections.Counter(["-","-","-","-","-","-","-"]) or collections.Counter(original_scheme[5]) == collections.Counter(["-","-","-","-","-","-","-","-"]) or collections.Counter(original_scheme[5]) == collections.Counter(["-"]):
					original_scheme[2] = "-"
				db_name = original_scheme[1]
				if "(" in db_name:
					new_db_name = db_name.split("(")[0]
				if lookup_allele_profile:
					#new_db_name=convert_mlst_to_pubMLST.convert(db_name)
					#unicode_name_for_lookup=bytestring.decode(new_db_name)
					#print("Downloading profile file for:", new_db_name)
					#profile_file = download_MLST_files(new_db_name, docFile)
					profile_file = mlst_db_path+"/pubmlst/"+db_name+"/"+db_name+".txt"
					print("Looking up:", current_alleles, "in", profile_file)
					db_file = open(profile_file,'r')
					db_line=db_file.readline().strip()
					while db_line is not None and db_line != '':
						db_alleles = db_line.split("\t")
						match = False
						db_type = db_alleles[0]
						for i in range(1,len(current_alleles)):
							#print(i, len(db_alleles), len(current_alleles))
							if db_alleles[i+1] == current_alleles[i]:
								match = True
							else:
								match = False
								break
						if match is False:
							db_line = db_file.readline().strip()
						else:
							print("matched:", db_line)
							original_scheme[2] = db_type
							break
				print("Expecting to add: ", original_scheme)
				print("B",checked_schemes)
				checked_schemes.append(original_scheme)
				print("A",checked_schemes)
			print("222-",original_scheme)
	for i in checked_schemes:
		print(i)

	checked_schemes.sort(reverse=True, key = lambda x: x[2])

	for i in checked_schemes:
		#print(i[1], i[5][2])
		if i[1] == "abaumannii" or i[1] == "abaumannii(Oxford)" or i[1] == "Acinetobacter_baumannii#1":
			if str(i[5][2]) == "182" or str(i[5][2]) == "189":
				print("FAKE NEWS!!!")
				i[2] = str(i[2])+"-PARALOG"
			i[1] = "abaumannii(Oxford)"
		elif i[1] == "abaumannii_2" or i[1] == "Acinetobacter_baumannii#2":
			i[1] = "abaumannii(Pasteur)"
		# I know this looks backwards....its not
		elif i[1] == "ecoli_2" or i[1] == "ecoli_achtman_4":
			i[1] = "ecoli(Achtman)"
		elif i[1] == "ecoli" or i[1] == "Escherichia_coli#2":
			i[1] = "ecoli_2(Pasteur)"


	# Print and check sets manually
	for i in checked_schemes:
		print("x",i)

	# Consolidate identical novel allele sets (This is not capable of doing more than 2 right now....but I dont see that eevr happening, so)
	checked_and_deduped_schemes=[]
	novel_allele_sets=[]
	for i in checked_schemes:
		if i[2] != "Novel_allele":
			checked_and_deduped_schemes.append(i)
		else:
			novel_allele_sets.append(i)

	for i in range(0,len(novel_allele_sets)):
		if len(novel_allele_sets) == 0:
			break
		if len(novel_allele_sets) == 1:
			checked_and_deduped_schemes.append(novel_allele_sets[0])
			novel_allele_sets.pop(0)
			break
		print(len(novel_allele_sets), novel_allele_sets)
		primary_db = novel_allele_sets[i][1]
		primary_alleles = []
		primary_source = novel_allele_sets[i][6]
		match_found=False
		for allele in novel_allele_sets[i][5]:
			if '*' in allele or '?' in allele or '~' in allele or '-' in allele:
				# Mark allele as found, but a mismatch
				primary_alleles.append('Mismatch')
			elif '-' in allele:
				# Mark allele as not found
				primary_alleles.append('Not Found')
			else:
				# Mark allele as found and complete
				primary_alleles.append('Complete')
		for j in range(0,len(novel_allele_sets)):
			if i == j:
				next
			else:
				if primary_db == novel_allele_sets[j][1]:
					secondary_alleles=[]
					secondary_source = novel_allele_sets[j][6]
					for allele in novel_allele_sets[j][5]:
						if '*' in allele or '?' in allele or '~' in allele or '-' in allele:
							# Mark allele as found, but a mismatch
							secondary_alleles.append('Mismatch')
						elif '-' in allele:
							# Mark allele as not found
							secondary_alleles.append('Not Found')
						else:
							# Mark allele as found and complete
							secondary_alleles.append('Complete')
					if primary_alleles == secondary_alleles:
						if primary_source == "standard":
							if secondary_source == "standard":
								print("Weird, shouldnt have both novel allele sets in one database come from assembly MLST")
							elif secondary_source == "standard/srst2":
									print("Weird, shouldnt have both novel allele sets in one database come from assembly MLST (plus a confirmation srst2)")
							elif secondary_source == "srst2":
								new_source = "standard/srst2"
							else:
								print("Weird, dont know what the secondary source is")
							primary_is_srst2=True
						elif primary_source == "srst2":
							if secondary_source == "standard":
								new_source = "standard/srst2"
							elif secondary_source == "standard/srst2":
								print("Weird, shouldnt have both novel allele sets in one database come from srst2 MLST (plus a confirmation assembly MLST)")
							elif secondary_source == "srst2":
								print("Weird, shouldnt have both novel allele sets in one database come from srst2")
							else:
								print("Weird, dont know what the secondary source is")
							secondary_is_srst2=True
						new_allele_list=[]
						for k in range(0,len(primary_alleles)):
							if primary_alleles[k] == 'Mismatch':
								if primary_is_srst2:
									new_allele_list.append('^'+novel_allele_sets[i][5][k]+','+novel_allele_sets[j][5][k])
								elif secondary_is_srst2:
									new_allele_list.append(novel_allele_sets[i][5][k]+',^'+novel_allele_sets[j][5][k])
								else:
									print("Weird, there is no srst2 found when comparing....these are all assembly MLSTs")
							else:
								new_allele_list.append(novel_allele_sets[j][5][k])
						new_entry = [novel_allele_sets[i][0], novel_allele_sets[i][1],novel_allele_sets[i][2],novel_allele_sets[i][3],novel_allele_sets[i][4],new_allele_list,new_source,novel_allele_sets[i][7]]
						checked_and_deduped_schemes.append(new_entry)
						novel_allele_sets.pop(j)
						novel_allele_sets.pop(i)
						match_found=True
				else:
					# Databases do not match, dont want to consolidate novel allele sets fron different DBs
					next
		if not match_found:
			checked_and_deduped_schemes.append(novel_allele_sets[i])
			novel_allele_sets.pop(i)

	all_Types_are_complete="unknown"
	outfile=isolate_name+"_combined.tsv"
	print(outfile)
	with open(outfile,'w') as writer:
		writer.write("Sample\tSource\tPulled on\tDatabase\tST\tlocus_1\tlocus_2\tlocus_3\tlocus_4\tlocus_5\tlocus_6\tlocus_7\tlocus_8\tlocus_9\tlocus_10\n")
		if len(checked_and_deduped_schemes) == 0:
			print("No schemes found")
			writer.write(isolate_name+"\tNone-"+genus+" "+species+"\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t")
			all_Types_are_complete="False"
		else:
			for i in checked_and_deduped_schemes:
				#print(i)
				allele_section=""
				if i[3] == 1:
					allele_section="-\n"
				else:
					for j in range (0,len(i[4])):
						if len(i[4][j].split("_")) > 1:
							gene_id=i[4][j].split("_")[1]
						else:
							gene_id=i[4][j]
						allele_section=allele_section+gene_id+"("+i[5][j]+")	"
					allele_section.strip()
					allele_section=allele_section+"\n"
				# Checks if current ST is either a number, has '-PARALOG' in it and that all previous ST types have completed
				if (i[2].isnumeric() or "-PARALOG" in i[2] or "Novel_profile" in i[2]) and i[2] != "-" and all_Types_are_complete != "False":
					#print(i[2], "= good")
					all_Types_are_complete = "True"
				else:
					#print(i[2], "= bad")
					all_Types_are_complete = "False"
				writer.write(isolate_name+"\t"+i[6]+"\t"+i[7]+"\t"+i[1]+"\t"+i[2]+"\t"+allele_section)
	status_file=isolate_name+"_status.txt"
	with open(status_file,'w') as writer:
		writer.write(all_Types_are_complete+"\n")



def main():
	print("Parsing MLST file ...")
	args = parseArgs()
	profile_lines=[]
	if args.input is not None:
		if os.stat(args.input).st_size > 0:
			#input_files=[[args.input , "standard"]]
			reg_file = open(args.input, 'r')
			### No headers complicate the difference between reg and srst2
			counter = 0
			for line in reg_file:
				print("reg:"+str(counter), line.replace("\n", ""))
				if counter > 0:
					print("appending -", line.replace("\n",""))
					profile_lines.append([line.replace("\n", ""), "mlst", args.input])
				counter+=1
		else:
			print("Input mlst file is empty")
	else:
		print("No regular mlst input file provided")
	if args.srst2 is not None:
		if os.stat(args.srst2).st_size > 0:
			srst2_file = open(args.srst2, 'r')
			counter = 0
			for line in srst2_file:
				print("srst2:"+str(counter), line.replace("\n", ""))
				if counter > 0:
					print("appending -", line.replace("\n",""))
					profile_lines.append([line.replace("\n", ""), "srst2", args.srst2])
				counter+=1
		else:
			print("Input SRST2 mlst file is empty")
	else:
		print("No srst2 input file provided")
	if len(profile_lines) > 0 and args.taxonomy is not None:
		do_MLST_check(profile_lines, args.taxonomy, args.mlst_database)
	else:
		print("No mlst files to check and fix")

if __name__ == '__main__':
	main()
