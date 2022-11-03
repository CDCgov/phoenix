#!/usr/bin/env python

#
# Description: Script to go through an isolates MLST output and check and reformat any instances of partial or multiple profiles found. SUB for profile or allele means that there is
# 	a novel allele or the profile has not been assigned yet, hence something needs to be submitted to pubmlst. AU (Allele Unknown) implies that an allele can not be determined
# 	no profile can be determined from the current assembly/qc_reads
#
# Usage: is python3 ./check_and_fix_MLST.py [-i input_torsten_MLST_file] [-s input_srst2_mlst_file] -t taxonomy_file
#
# Modules required: None
#
# v2.0 (9/23/2022)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#from asyncio.windows_events import NULL
from locale import currency
import sys
import os
import glob
import math
import copy
import argparse
import collections
import itertools as it
#from pathlib import Path
import convert_mlst_to_pubMLST
import xml.dom.minidom as xml
import urllib2 as url
import re
from urlparse import urlparse
import ssl
from datetime import datetime
from os import path


db_lookup = {}

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to check MLST types for duplicate alleles and implications on final typing')
	parser.add_argument('-i', '--input', required=False, help='input mlst filename')
	parser.add_argument('-s', '--srst2', required=False, help='srst2 input file')
	parser.add_argument('-t', '--taxonomy', required=True, help='Location of taxonomy file to pull right scheme')
	return parser.parse_args()

# test if a node is an Element and that it has a specific tag name
def testElementTag(node, name):
		return node.nodeType == node.ELEMENT_NODE and node.localName == name

# Get the text from an element node with a text node child
def getText(element):
	result = ''
	for node in element.childNodes:
		if node.nodeType == node.TEXT_NODE:
			result += node.data
	return normaliseText(result)

# remove unwanted whitespace including linebreaks etc.
def normaliseText(str):
	return ' '.join(str.split())

# A collection of interesting information about a taxa
class SpeciesInfo(object):
	def __init__(self):
		self.name = None # String name of species
		self.database_url = None # URL as string
		self.retrieved = None # date as string
		self.profiles_url = None # URL as string
		self.profiles_count = None # positive integer
		self.loci = [] # list of loci

class LocusInfo(object):
	def __init__(self):
		self.url = None
		self.name = None

# retrieve the interesting information for a given sample element
def getSpeciesInfo(species_node, species, exact):
	this_name2 = getText(species_node)
	this_name = getText(species_node).encode('utf-8')
	store = False
	#print(this_name, species, exact)
	#print(this_name2, species, exact)

	# if this_name == species:
	# 	print "Match"
	# else:
	# 	print "No match"
	#
	if exact:
		if this_name == species:
			store = True
	else:
		if this_name.startswith(species):
			store = True
		#else:
		#	print "says not True with starts with"
	if store:
	#if True:
		info = SpeciesInfo()
		info.name = this_name
		for mlst_node in species_node.getElementsByTagName('mlst'):
			for database_node in mlst_node.getElementsByTagName('database'):
				for database_child_node in database_node.childNodes:
					if testElementTag(database_child_node, 'url'):
						info.database_url = getText(database_child_node).encode('utf-8')
					elif testElementTag(database_child_node, 'retrieved'):
						info.retrieved = getText(database_child_node).encode('utf-8')
					elif testElementTag(database_child_node, 'profiles'):
						for profile_count in database_child_node.getElementsByTagName('count'):
							info.profiles_count = getText(profile_count).encode('utf-8')
						for profile_url in database_child_node.getElementsByTagName('url'):
							info.profiles_url = getText(profile_url).encode('utf-8')
					elif testElementTag(database_child_node, 'loci'):
						for locus_node in database_child_node.getElementsByTagName('locus'):
							locus_info = LocusInfo()
							locus_info.name = getText(locus_node).encode('utf-8')
							for locus_url in locus_node.getElementsByTagName('url'):
								locus_info.url = getText(locus_url).encode('utf-8')
							info.loci.append(locus_info)
		return info
	else:
		return None

def download_MLST_files(tax_to_download):
	ssl._create_default_https_context = ssl._create_unverified_context
	docFile = url.urlopen("http://pubmlst.org/data/dbases.xml")
	force=False
	if tax_to_download == "Streptococus thermophilus":
		force=True
	doc = xml.parse(docFile)
	root = doc.childNodes[0]
	found_species = []
	for species_node in root.getElementsByTagName('species'):
		info = getSpeciesInfo(species_node, tax_to_download, force)
		if info != None:
			found_species.append(info)
	if len(found_species) == 0:
		print ("No species matched your query.")
		exit(1)

	#print(len(found_species))
	assert len(found_species) == 1
	species_info = found_species[0]
	species_name_underscores = species_info.name.replace(' ', '_')
	species_name_underscores = species_name_underscores.replace('/', '_')


	# output information for the single matching species
	#species_all_fasta_filename = work_dir+species_name_underscores + '.fasta'
	#species_all_fasta_file = open(species_all_fasta_filename, 'w')
	#print(type(work_dir), work_dir)
	log_filename = "mlst_data_download_{}_{}.log".format(species_name_underscores, species_info.retrieved)
	log_file = open(log_filename, "w")
	profile_path = urlparse(species_info.profiles_url).path.encode('utf-8')
	#print(type(work_dir), work_dir)
	#print(type(profile_path), profile_path)
	profile_filename = species_name_underscores+"_"+profile_path.split('/')[-1].replace("profiles_csv", "profiles.csv")
	#log_file.write("definitions: {}\n".format(profile_filename))
	log_file.write("{} profiles\n".format(species_info.profiles_count))
	#log_file.write("sourced from: {}\n\n".format(species_info.profiles_url))
	profile_doc = url.urlopen(species_info.profiles_url)
	profile_file = open(profile_filename, 'w')
	profile_file.write(profile_doc.read())
	profile_file.close()
	profile_doc.close()
	# for locus in species_info.loci:
	# 	locus_path = urlparse(locus.url).path
	# 	locus_filename = work_dir+locus_path.split('/')[-1]
	# 	log_file.write("locus {}\n".format(locus.name))
	# 	log_file.write(locus_filename + '\n')
	# 	log_file.write("Sourced from {}\n\n".format(locus.url))
	# 	locus_doc = url.urlopen(locus.url)
	# 	locus_file = open(locus_filename, 'w')
	# 	locus_fasta_content = locus_doc.read()
	# 	locus_file.write(locus_fasta_content)
	# 	species_all_fasta_file.write(locus_fasta_content)
	# 	locus_file.close()
	# 	locus_doc.close()
	# log_file.write("all loci: {}\n".format(species_all_fasta_filename))
	log_file.close()
	#species_all_fasta_file.close()

	# print "\n  For SRST2, remember to check what separator is being used in this allele database"
	# head = os.popen('head -n 1 ' + species_all_fasta_filename).read().rstrip()
	# m = re.match('>(.*)([_-])(\d*)',head).groups()
	# if len(m)==3:
	# 	print
	# 	print "  Looks like --mlst_delimiter '" + m[1] + "'"
	# 	print
	# 	print "  " + head + "  --> -->  ",
	# 	print m
	# print
	# print "  Suggested srst2 command for use with this MLST database:"
	# print
	# print "    srst2 --output test --input_pe *.fastq.gz --mlst_db " + species_name_underscores + '.fasta',
	# print "--mlst_definitions " + format(profile_filename),
	# print "--mlst_delimiter '" + m[1] + "'"
	today=datetime.today().strftime('%Y-%m-%d')
	if path.isfile(species_name_underscores+"_pull_dates.txt"):
		with open(species_name_underscores+"_pull_dates.txt",'a+') as pull_file:
			dates=pull_file.read()
			if today in dates:
				print("Today is already listed in pull file, not appending")
			else:
				print("Today is not in the list of pull dates, appending today")
				pull_file.write(today+"\n")
	else:
		with open(species_name_underscores+"_pull_dates.txt",'w') as pull_file:
			print("pull dates file does not exist, creating and appending today")
			pull_file.write(today)



	return profile_filename

# main function that looks if all MLST types are defined for an outptu mlst file
def do_MLST_check(input_MLST_line_tuples, taxonomy_file):
	location="/".join(taxonomy_file.split("/")[0:-1])+"/mlst"
	#print taxonomy_file
	isolate_name = taxonomy_file.split(".")[:-1]
	isolate_name = ".".join(isolate_name)
	tax_file = open(taxonomy_file, 'r')
	today=datetime.today().strftime('%Y-%m-%d')
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
			original_items=original_line.split("	")
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
	#					print(":", original_items[allele])
						allele_Identifier=original_items[allele].split("(")[0]
						alleles=original_items[allele].split("(")[1].split(")")[0].split(",")
						if len(alleles) > 1:
							original_type="-"
						allele_names.append(allele_Identifier)
						allele_list.append(alleles)

					# Create a template to build Scheme array
					expanded_allele_list.append([sample,db_name, original_type, len(allele_names), allele_names, [], "standard", today])

					# Test first parse
	#				print("e:",expanded_allele_list)
	#				print("o:", allele_list)

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
	#								print("0 - Going to add:", allele[i])
									for j in expanded_allele_list:
	#									print("to - ", j)
										j[5].append(allele[i])
	#									print("After adding - ", j)
								else:
									temp2_allele_list=copy.deepcopy(temp_allele_list)
	#								print(temp2_allele_list)
	#								print(i,"- Going to add:", allele[i])
									for j in temp2_allele_list:
	#									print("to - ", j)
										j[5].append(allele[i])
	#									print("After adding - ", j)
									expanded_allele_list.append(copy.deepcopy(temp2_allele_list)[0])
	#					print(expanded_allele_list)
				else:
					#expanded_allele_list.append([sample,db_name, "-", 7, ["-","-","-","-","-","-","-"], ["-","-","-","-","-","-","-"], "standard"])
					expanded_allele_list.append([sample,db_name, "-", 1, ["-"], ["-"], "standard", today])


			#allele_list=[['1'], ['3'], ['189','3'], ['2'], ['2'], ['96','107'], ['3']]
			elif MLST_filetype == "srst2":
				print original_items,len(original_items)
				sample=original_items[0]
				pubmlst_db_name=original_items[1]
				db_name=convert_mlst_to_pubMLST.back_2_MLST(pubmlst_db_name)
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

					expanded_allele_list.append([sample,db_name, original_type, len(allele_names), allele_names, [], "srst2", today])

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
	#								print("0 - Going to add:", allele[i])
									for j in expanded_allele_list:
	#									print("to - ", j)
										j[5].append(allele[i])
	#									print("After adding - ", j)
								else:
									temp2_allele_list=copy.deepcopy(temp_allele_list)
	#								print(temp2_allele_list)
	#								print(i,"- Going to add:", allele[i])
									for j in temp2_allele_list:
	#									print("to - ", j)
										j[5].append(allele[i])
	#									print("After adding - ", j)
									expanded_allele_list.append(copy.deepcopy(temp2_allele_list)[0])
						print(expanded_allele_list)
				else:
					expanded_allele_list.append([sample,db_name, "-", 1, ["-"], ["-"], "srst2", today])

			else:
				print("Unknown MLST filetype, can not continue")
				exit()
			#original_line=MLST_file.readline().strip("	")
			original_schemes.append(expanded_allele_list)

		#MLST_file.close()

	### Shows the list of what schemes were found
	print("Schemes found:", len(original_schemes))
	catted_scheme_list=[]
	for oscheme in original_schemes:
		print oscheme
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
	#		print("Checking:", catted_scheme_list[i][5], catted_scheme_list[j][5])
			#print catted_scheme_list[i],catted_scheme_list[j]
			if collections.Counter(catted_scheme_list[i][5]) == collections.Counter(catted_scheme_list[j][5]):
				print("SAME!!!!", i, j)
				if catted_scheme_list[min(i,j)][6] != catted_scheme_list[max(i,j)][6]:
					print "a"
					if catted_scheme_list[min(i,j)][6] == "standard":
						print "b"
						if catted_scheme_list[max(i,j)][6] == "srst2":
							print "c"
							#catted_scheme_list[min(i,j)][4] = catted_scheme_list[max(i,j)][4]
							catted_scheme_list[min(i,j)][6] = "standard/srst2"
							dupes.append(max(i,j))
						elif catted_scheme_list[max(i,j)][6] == "standard/srst2":
							print "d"
							dupes.append(min(i,j))
					elif catted_scheme_list[min(i,j)][6] == "srst2":
						print "e"
						if catted_scheme_list[max(i,j)][6] == "standard":
							print "f"
							#catted_scheme_list[max(i,j)][4] = catted_scheme_list[min(i,j)][4]
							catted_scheme_list[max(i,j)][6] = "standard/srst2"
							dupes.append(min(i,j))
						elif catted_scheme_list[max(i,j)][6] == "standard/srst2":
							print "g"
							dupes.append(min(i,j))
					elif catted_scheme_list[min(i,j)][6] == "standard/srst2":
						print "h"
						dupes.append(max(i,j))
					else:
						print "i"
						print("Should never have something that is not standard, srst2, or standard/srst2")
				else:
					if max(i,j) not in dupes:
						print "j"
						dupes.append(max(i,j))
					else:
						print("Index already Found")
				print(dupes)
			else:
	#			print(catted_scheme_list[i][5], "does not equal", catted_scheme_list[i][5])
				continue
	print "k"
	print(dupes)
	dedupped_dupes=list(set(dupes))
	dedupped_dupes.sort(reverse=True)
	print(dedupped_dupes)
	print(catted_scheme_list)
	for k in dedupped_dupes:
		print k
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
					#	MLST
					#	~ : full length novel allele
					#	? : partial match (>min_cov & > min_ID). Default min_cov = 10, Default min_ID=95%
					#	- : Allele is missing
					#
					#	srst2
					#	* : Full length match with 1+ SNP. Novel
					#	? : edge depth is below N or average depth is below X. Default edge_depth = 2, Default average_depth = 5
					#	- : No allele assigned, usually because no alleles achieved >90% coverage



					if '*' in allele or '?' in allele or '~' in allele or '-' in allele:
						original_scheme[2] = "1+_Novel_allele"
						bad_alleles += 1
						lookup_allele_profile = False
				## Reapply dash showing no db or proximal scheme was ever found
				if collections.Counter(original_scheme[5]) == collections.Counter(["-","-","-","-","-","-","-"]) or collections.Counter(original_scheme[5]) == collections.Counter(["-","-","-","-","-","-","-","-"]) or collections.Counter(original_scheme[5]) == collections.Counter(["-"]):
					original_scheme[2] = "-"
				db_name = original_scheme[1]
				if "(" in db_name:
					new_db_name = db_name.split("(")[0]
				if lookup_allele_profile:
					new_db_name=convert_mlst_to_pubMLST.convert(db_name)
					#unicode_name_for_lookup=bytestring.decode(new_db_name)
					print("Downloading profile file for:", new_db_name)
					profile_file = download_MLST_files(new_db_name)
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


	for i in checked_schemes:
		print(i)

	outfile=isolate_name+"_combined.tsv"
	print(outfile)
	with open(outfile,'w') as writer:
		writer.write("Sample	Source	Pulled on	Database	ST	locus_1	locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 locus_8 locus_9 locus_10\n")
		if len(checked_schemes) == 0:
			print("No schemes found")
			writer.write(isolate_name+"	None-"+genus+" "+species+"	-")
		else:
			for i in checked_schemes:
				print(i)
				allele_section=""
				if i[3] == 1:
					allele_section="-"
				else:
					for j in range (0,len(i[4])):
						if len(i[4][j].split("_")) > 1:
							gene_id=i[4][j].split("_")[1]
						else:
							gene_id=i[4][j]
						allele_section=allele_section+gene_id+"("+i[5][j]+")	"
					allele_section.strip()
					allele_section=allele_section+"\n"
				writer.write(isolate_name+"	"+i[6]+"	"+i[7]+"	"+i[1]+"	"+i[2]+"	"+allele_section)



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
	print("No srst2 input file provided")
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
	do_MLST_check(profile_lines, args.taxonomy)
else:
	print("No mlst files to check and fix")
