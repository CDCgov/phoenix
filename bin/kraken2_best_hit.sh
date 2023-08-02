#!/bin/bash -l

#
# Description: Grabs the best species match based on %/read hits from the kraken tool run. Simplified for nextflow inclusion
#
# Usage: ./kraken_best_hit.sh -i path_to_.list_file
#
# Output location: same as input path_to_.list_file
#
# Modules required: None
#
# v1.0.3 (04/19/2022)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./kraken_best_hit.sh -i path_to_list_file [-q count_file (reads or congtigs)] -n sample_name"
	echo "Output is saved to folder where .list file exists"
}

# Parse command line options
options_found=0
while getopts ":h?i:q:n:t:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		i)
			echo "Option -i triggered, argument = ${OPTARG}"
			list_file=${OPTARG}
			;;
		q)
			echo "Option -q triggered, argument = ${OPTARG}"
			read_file=${OPTARG}
			;;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			sample_name=${OPTARG}
			;;
		t)
			echo "Option -t triggered"
			terra=${OPTARG}
			;;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

if [[ "${options_found}" -eq 0 ]]; then
	echo "No argument supplied to kraken2_best_hit.sh, exiting"
	show_help
	exit 1
fi

# set the correct path for bc - needed for terra
if [[ $terra = "terra" ]]; then
	bc_path=/opt/conda/envs/phoenix/bin/bc
else
	bc_path=bc
fi

# Based upon standard naming protocols pulling last portion of path off should result in proper name
#sample_name=$(echo "${list_file}" | rev | cut -d'/' -f1 | cut -d'.' -f4- | rev)
#SAMPDATADIR=$(echo "${list_file}" | rev | cut -d'/' -f2- | rev)
#sample_name=$(basename "${list_file}" .kraken2.summary.txt)
#echo "-${sample_name}-"

#Creates the default values for output in case any calculations are interrupted.
unclass_reads=0
class_reads=0
u_percent=0
unclass_percent=0
domain_percent=0
phylum_percent=0
class_percent=0
order_percent=0
family_percent=0
genus_percent=0
species_percent=0
domain_percent_total=0
phylum_percent_total=0
class_percent_total=0
order_percent_total=0
family_percent_total=0
genus_percent_total=0
species_percent_total=0
domain_reads=0
phylum_reads=0
class_reads=0
order_reads=0
family_reads=0
genus_reads=0
species_reads=0
classified_reads=0
domain="N/A"
phylum="N/A"
class="N/A"
order="N/A"
family="N/A"
genus="N/A"
species="N/A"

#Parses the kraken output list line by line
while IFS= read -r line  || [ -n "$line" ]; do
	# Removes the whitespace at the end of the line
	line=${line}
	# Turns the line into an array delimited by whitespace
	arrLine=(${line})
	# First element in array is the percent of reads identified as the current taxa
	percent=${arrLine[0]}
	# 3rd element is the taxon level classification
	classification=${arrLine[3]}
	# 2nd element is the actual read count identified as the current taxa
	reads=${arrLine[1]}
	# All elements beyond 5 are stored as the description of the taxa, genus/species/strain and any other included information
	description="${arrLine[*]:5}"
	# Pulls the first word from the description to find if it is root, the descriptor used to identify any number of reads classified, but higher than domain
	first_desc=$(echo "${description}" | cut -d' ' -f1)
	# Unclassified read info (reads and percent) is captured if the classifier claims to be 'u'
	# Kraken lists entries in order of abundance (highest first) therefore each time a new taxon level is encountered it will be the highest ranking representative

	#echo "${line}---${percent}---${classification}---${reads}---${description}"

	#echo "${reads} ${domain_reads}"
	if [ "${classification}" = "U" ]; then
		unclass_percent=${percent}
		# Will be overwritten by pre flag, but not in post
		unclass_reads=${reads}
		# Grabs all supra domain level read info (reads and percent)
	elif [ "${classification}" = "R" ] || ([ "${classification}" = "-" ] && [ "${first_desc}" = "root" ]); then
		classified_reads="${reads}"
		root_percent="${percent}"
	# Grabs all read info (identifier, reads and percent) for best domain level entry
	elif [ "${classification}" = "D" ] && [ "${reads}" -gt "${domain_reads}" ]; then
		domain=${description^}
		domain_percent=${percent}
		domain_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best phylum level entry
	elif [ "${classification}" = "P" ] && [ "${reads}" -gt "${phylum_reads}" ]; then
		phylum=${description^}
		phylum_percent=${percent}
		phylum_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best class level entry
	elif [ "${classification}" = "C" ] && [ "${reads}" -gt "${class_reads}" ]; then
		class=${description^}
		class_percent=${percent}
		class_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best order level entry
	elif [ "${classification}" = "O" ] && [ "${reads}" -gt "${order_reads}" ]; then
		order=${description^}
		order_percent=${percent}
		order_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best family level entry
	elif [ "${classification}" = "F" ] && [ "${reads}" -gt "${family_reads}" ]; then
		family=${description^}
		family_percent=${percent}
		family_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best genus level entry
	elif [ "${classification}" = "G" ] && [ "${reads}" -gt "${genus_reads}" ]; then
		genus=${description^}
		genus_percent=${percent}
		genus_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best species level entry
	elif [ "${classification}" = "S" ] && [ "${reads}" -gt "${species_reads}" ]; then
		echo "Old: ${species}-${species_reads}"
		gs=(${description^})
		species=${gs[@]:1}
		species_percent=${percent}
		species_reads=${reads}
		echo "New: ${species}-${species_reads}"
	fi
done < "${list_file}"

#echo -e "Test: ${unclass_reads}-${unclass_percent}\n${classified_reads}-${root_percent}\n${domain_reads}-${domain_percent}-${domain_percent_total}\n${domain_reads}-${domain_percent}-${domain_percent_total}\n${phylum_reads}-${phylum_percent}-${phylum_percent_total}\n${class_reads}-${class_percent}-${class_percent_total}\n${order_reads}-${order_percent}-${order_percent_total}\n${family_reads}-${family_percent}-${family_percent_total}\n${genus_reads}-${genus_percent}-${genus_percent_total}\n${species_reads}-${species_percent}-${species_percent_total}\n"

# ############ To BE determined as to exact naming conventions
# if [[ "${list_file}" = *"/kraken2_asmbld_weighted/"* ]]; then
# #	echo "${unclass_percent}:${root_percent}:${domain_percent}:${phylum_percent}:${class_percent}:${order_percent}:${family_percent}:${genus_percent}:${species_percent}"
# 	total_percent=$(echo "${unclass_percent} + ${root_percent}" | bc)
# 	unclass_percent=$(echo "${unclass_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# 	domain_percent=$(echo "${domain_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# 	phylum_percent=$(echo "${phylum_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# 	order_percent=$(echo "${order_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# 	class_percent=$(echo "${class_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# 	family_percent=$(echo "${family_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# 	genus_percent=$(echo "${genus_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# 	species_percent=$(echo "${species_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# #	echo "${unclass_percent}:${root_percent}:${domain_percent}:${phylum_percent}:${class_percent}:${order_percent}:${family_percent}:${genus_percent}:${species_percent}"
# fi

############################### NOT CURRENTLY IN USE ##########################################################################################

# Calculate % of unclassified reads/contigs using sum of highest taxon level reads against total reads found in QC counts
# Grabs total possible reads from preQC counts if kraken was used on reads (pre assembly)
if [[ "${list_file}" = *"kraken2_trimd.summary.txt" ]]; then
	echo "doing trimd"
#	r1s=$(tail -n 1 "${read_file}" | cut -d'	' -f2)
#	r2s=$(tail -n 1 "${read_file}" | cut -d'	' -f4)
#	file_reads=$(( r1s + r2s ))
	# Calculates the true count of unclassified reads/contigs rather than the reported value from kraken
#	unclass_reads=$(( file_reads - classified_reads ))
	# Calculates the percent of unclassified reads/contigs using the total possible reads
#	u_percent=$(echo "${unclass_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# Grabs total possible bases from contigs in trimmed assembly (post assembly, using weighted kraken output)
elif [[ "${list_file}" = *"kraken2_asmbld.summary.txt" ]]; then
	echo "doing asmbld"
	# Full length of assembly? Still not 100% sure how the kreport uses read lengths
#	file_reads=$(head -n 14 "${read_file}" | tail -n1 | cut -d$'\t' -f2)
	# Calculates percent of classified reads as 100*classified reads/contigs
#	u_percent=$(echo "${unclass_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
elif [[ "${list_file}" = *".kraken2_wtasmbld.summary.txt" ]]; then
	echo "doing weighted"
	# Full length of assembly? Still not 100% sure how the kreport uses read lengths
#	file_reads=$(head -n 14 "${read_file}" | tail -n1 | cut -d$'\t' -f2)
	# Calculates percent of classified reads as 100*classified reads/contigs
#	u_percent=$(echo "${unclass_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	#echo "${unclass_percent}:${root_percent}:${domain_percent}:${phylum_percent}:${class_percent}:${order_percent}:${family_percent}:${genus_percent}:${species_percent}"
	total_percent=$(echo "${unclass_percent} + ${root_percent}" | $bc_path)
	unclass_percent=$(echo "${unclass_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	root_percent=$(echo "${root_percent} ${total_percent}" | awk '{ printf "%2.2f" , ($1*100)/$2 }' )
	domain_percent=$(echo "${domain_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	phylum_percent=$(echo "${phylum_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	order_percent=$(echo "${order_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	class_percent=$(echo "${class_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	family_percent=$(echo "${family_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	genus_percent=$(echo "${genus_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	species_percent=$(echo "${species_percent} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	#echo "${unclass_percent}:${root_percent}:${domain_percent}:${phylum_percent}:${class_percent}:${order_percent}:${family_percent}:${genus_percent}:${species_percent}"
fi


#Print out the best taxa for each level and its corresponding % of reads reported by kraken, % reads of total, taxon description
# echo -e "U: ${unclass_percent} unclassified\\nD: ${domain_percent} ${domain}\\nP: ${phylum_percent} ${phylum}\\nC: ${class_percent} ${class}\\nO: ${order_percent} ${order}\\nF: ${family_percent} ${family}\\nG: ${genus_percent} ${genus}\\ns: ${species_percent} ${species}" > "${sample_name}.summary.txt"

###With headers
echo -e "Taxon level	Match percentage	Taxa\nU: ${unclass_percent} unclassified\\nD: ${domain_percent} ${domain}\\nP: ${phylum_percent} ${phylum}\\nC: ${class_percent} ${class}\\nO: ${order_percent} ${order}\\nF: ${family_percent} ${family}\\nG: ${genus_percent} ${genus}\\ns: ${species_percent} ${species}" > "${sample_name}.summary.txt"


#Script exited gracefully (unless something else inside failed)
exit 0
