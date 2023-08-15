#!/bin/bash -l
#
# Description: Creates a single file that attempts to pull the best taxonomic information from the isolate. Currently, it operates in a linear fashion, e.g. 1.ANI, 2.kraken2
# 	The taxon is chosen based on the highest ranked classifier first
#
# Usage: ./determine_texID.sh -k weighted_kraken_report -s sample_name -f formatted_fastani_file -d nodes_file -m names_file
#
# Modules required: None
#
# v2.0 (08/15/2023)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./determine_taxID.sh -k weighted_kraken_report -s sample_name -f formatted_fastani_file -d nodes_X.dmp -m names_X.dmp"
	echo "Output is saved to /sample_name/sample_name.tax"
}

# Parse command line options
options_found=0
while getopts ":h?k:s:f:d:r:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
			show_help
			exit 0
			;;
		k)
			echo "Option -k triggered, argument = ${OPTARG}"
			weighted_kraken=${OPTARG};;
		r)
			echo "Option -r triggered, argument = ${OPTARG}"
			trimmed_kraken=${OPTARG};;
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
		f)
			echo "Option -f triggered, argument = ${OPTARG}"
			fastani_file=${OPTARG};;
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			nodes=${OPTARG};;
   		m)
			echo "Option -m triggered, argument = ${OPTARG}"
			names=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

# Set default values for all taxonomic levels
Domain="Not_assigned"
Phylum="Not_assigned"
Class="Not_assigned"
Order="Not_assigned"
Family="Not_assigned"
Genus="Not_assigned"
species="Not_assigned"
source="Not_assigned"
confidence_index="0"
source_file="Not_assigned"


# Function to check which source to use as the 'determinator'. Single int parameter can be used to tell which level to jump in at
Check_source() {
	start_at="${1}"
	if [[ "${start_at}" -le 1 ]]; then
#		for f in ${OUTDATADIR}/ANI/*; do
		if [[ "${fastani_file}" = *".fastANI"*".txt" ]]; then
			header=$(head -n1 ${fastani_file})
			if [[ ${header} != "No matching ANI database found for"* ]] && [[ ${header} != "No MASH hit found" ]] && [[ ${header} != "0.00%"* ]] ; then
				do_ANI
			return
			else
				Check_source 2
			fi
		fi
#		done
	fi
	if [[ "${start_at}" -le 2 ]]; then
#		if [[ -s "${OUTDATADIR}/kraken2_weighted/${sample_name}.summary.txt" ]]; then
		if [[ -s "${weighted_kraken}" ]]; then
			do_kraken2_assembly
		return
		fi
	fi
	if [[ "${start_at}" -le 3 ]]; then
#		if [[ -s "${OUTDATADIR}/kraken2_weighted/${sample_name}.summary.txt" ]]; then
		if [[ -s "${trimmed_kraken}" ]]; then
			do_kraken2_reads
		return
		fi
	fi
	echo "No ACCEPTABLE source found to determine taxonomy"
}

# Function to pull info from ANI output
do_ANI() {
	source="ANI_REFSEQ"
#	source_file=$(find "${OUTDATADIR}/ANI" -type f -name "${sample_name}.fastANI_*.txt" | sort -k3,3 -Vrt '.' | head -n1)
	source_file="${fastani_file}"
	if [[ -s "${fastani_file}" ]]; then
 		header=$(head -n 1 "${fastani_file}")
   		info=$(tail -n 1 "${fastani_file}")
		Genus=$(echo "${info}" | cut -d'	' -f3 | cut -d' ' -f1)
		species=$(echo "${info}" | cut -d'	' -f3 | cut -d' ' -f2- | sed 's/[][]//g')
		confidence_index=$(echo "${info}" | cut -d'	' -f1)
	else
		echo "source file (${fastani_file}) is empty"
	fi
}

# Function to pull info from kraken2 output based on assembly
do_kraken2_assembly() {
	source="kraken2_wtasmbld"
#	source_file="${OUTDATADIR}/kraken2_weighted/${sample_name}.summary.txt"
	source_file="${weighted_kraken}"
	#echo "${source}"
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		echo $line
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken2, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{print $3}')
		elif [ "${first}" = "G" ]
		then
			Genus=$(echo "${line}" | awk -F ' ' '{print $3}')
		fi
	done < "${source_file}"
	confidence_index=$(tail -n1 "${source_file}" | cut -d' ' -f2)
	#confidence_index="${confidence_index}"
}

# Function to pull info from kraken2 output based on assembly
do_kraken2_reads() {
	source="kraken2_trimmed"
#	source_file="${OUTDATADIR}/kraken2_weighted/${sample_name}.summary.txt"
	source_file="${trimmed_kraken}"
	#echo "${source}"
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken2, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{ print $3 }')
		elif [ "${first}" = "G" ]
		then
			Genus=$(echo "${line}" | awk -F ' ' '{print $3}')
		fi
	done < "${source_file}"
	confidence_index=$(tail -n1 "${source_file}" | cut -d' ' -f2)
	#confidence_index="${confidence_index}"
}

# Start the program by checking ALL sources
Check_source 0


# Check if species was assigned and get starting taxID
if [[ -n ${species} ]]; then
	species=$(echo ${species} | tr -d [:space:])
	for name_line in $(zgrep $"	|	${species}	|	" ${names}); do
 		taxID=$(echo "${name_line}" | cut -d$'\t' -f1)
		name=$(echo "${name_line}" | cut -d$'\t' -f3)
		unique_name=$(echo "${name_line}" | cut -d$'\t' -f5)
		name_class=$(echo "${name_line}" | cut -d$'\t' -f7)
		if [[ "${name_class}" = "scientific name" ]]; then
			species_taxID="${taxID}"
		fi
	done
 fi

 # See if we can at least start at genus level to fill in upper taxonomy
 if [[ -z "${species_taxID}" ]]; then
 	# Check if genus was assigned
 	Genus=$(echo ${Genus} | tr -d [:space:] | tr -d "[]")
 	if [[ -n "${Genus}" ]]; then
  		for name_line in $(zgrep $"	|	${species}	|	" ${names}); do
  			taxID=$(echo "${name_line}" | cut -d$'\t' -f1)
			name=$(echo "${name_line}" | cut -d$'\t' -f3)
			unique_name=$(echo "${name_line}" | cut -d$'\t' -f5)
			name_class=$(echo "${name_line}" | cut -d$'\t' -f7)
			if [[ "${name_class}" = "scientific name" ]]; then
				genus_taxID="${taxID}"
			fi
  		done
    	fi
fi

# Check if genus was assigned as peptoclostridium and relabel it as Clostridium for downstream analyses relying on this older naming convention
#if [[ ${Genus} == "Peptoclostridium" ]]; then
#	Genus="Clostridium"
#fi

#taxa_indices=( "species" "genus" "family" "order" "class" "phylum" "kingdom")
taxa_indices=( "kingdom" "phylum" "class" "order" "family" "genus" "species")

declare -A taxID_list=( [kingdom]="NA" [phylum]="NA" [class]="NA" [order]="NA" [family]="NA" [genus]="NA" [species]="NA")
declare -A tax_name_list=( [kingdom]="NA" [phylum]="NA" [class]="NA" [order]="NA" [family]="NA" [genus]="NA" [species]="NA")


if [[ -z "${species_taxID}" ]]; then
	if [[ -z "${genus_taxID}" ]]; then
		echo -e "${source}	${confidence_index}	${source_file}\nK:	Unknown\nP:	Unknown\nC:	Unknown\nO:	Unknown\nF:	Unknown\nG:	Unknown\ns:	Unknown\n" > "${sample_name}.tax"
		exit
	else
		max_counter=6
		taxID_list[species]=0
		tax_name_list[species]="Unknown"
  		taxID="${genus_taxID}"
	fi
else
 	max_counter=7
	taxID="${species_taxID}"
fi
 
counter=0

while [[ ${counter} -lt "${max_counter}" ]]; do
	index=$(( max_counter - 1 - counter ))
	node_line=$(zgrep "^${taxID}	|	" ${nodes})
	parent=$(echo "${node_line}" | cut -d$'\t' -f3)
	rank=$(echo "${node_line}" | cut -d$'\t' -f5)
	#echo "${counter},${taxID},${parent},${rank},${taxa_indices[${index}]}"
	if [[ "${rank}" = "${taxa_indices[${index}]}" ]] || [[ "${rank}" = "superkingdom" ]]; then
		#echo "Adding-${taxa_indices[${index}]}-${taxID}"
		taxID_list["${taxa_indices[${index}]}"]="${taxID}"
		counter=$(( counter + 1 ))
	fi
	taxID=${parent}
	#for taxa in "${!taxID_list[@]}"; do
	#	echo "${taxa}-${taxID_list[${taxa}]}"
	#done
done

counter=0
while [[ ${counter} -lt "${max_counter}" ]]; do
	#index=$(( max_counter - 1 - counter ))
	taxID=${taxID_list[${taxa_indices[${counter}]}]}
	IFS=$'\n'
	for name_line in $(zgrep "^${taxID}	|	" ${names}); do
		name=$(echo "${name_line}" | cut -d$'\t' -f3)
		unique_name=$(echo "${name_line}" | cut -d$'\t' -f5)
		name_class=$(echo "${name_line}" | cut -d$'\t' -f7)
		if [[ "${name_class}" = "scientific name" ]]; then
			#echo "${name^} as ${taxa_indices[${counter}]}"
			tax_name_list[${taxa_indices[${counter}]}]="${name^}"
		fi
	done
	counter=$(( counter + 1 ))
done

if [[ "${tax_name_list[species]}" != "Unknown" ]]; then
	tax_name_list[species]=$(echo ${tax_name_list[species],,} | cut -d' ' -f2-)
fi

for i in "${!tax_name_list[@]}"; do
	echo $i-${tax_name_list[$i]}
done

# Print output to tax file for sample
#echo -e "(${source})-${confidence_index}-${source_file}\nD:	${Domain}\nP:	${Phylum}\nC:	${Class}\nO:	${Order}\nF:	${Family}\nG:	${Genus}\ns:	${species}\n" > "${OUTDATADIR}/${sample_name}.tax"
echo -e "${source}	${confidence_index}	${source_file}\nK:${taxID_list[0]}	${taxID_list[0]}\nP:${taxID_list[1]}	${taxID_list[1]}\nC:${taxID_list[2]}	${taxID_list[2]}\nO:${taxID_list[3]}	${taxID_list[3]}\nF:${taxID_list[4]}	${taxID_list[4]}\nG:${taxID_list[5]}	${taxID_list[5]}\ns:${taxID_list[6]}	${taxID_list[06}\n" > "${sample_name}.tax"
