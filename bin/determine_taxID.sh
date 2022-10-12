#!/bin/bash -l
#
# Description: Creates a single file that attempts to pull the best taxonomic information from the isolate. Currently, it operates in a linear fashion, e.g. 1.ANI, 2.kraken2
# 	The taxon is chosen based on the highest ranked classifier first
#
# Usage: ./determine_texID.sh -k weighted_kraken_report -s sample_name -f formatted_fastani_file -d database_file(taxes.csv location)
#
# Modules required: None
#
# v1.0.8 (08/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./determine_taxID.sh -k weighted_kraken_report -s sample_name -f formatted_fastani_file -d database_file(taxes.csv location)"
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
			tax_DB=${OPTARG};;
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
			if [[ ${header} != "No matching ANI database found for"* ]] && [[ ${header} != "0.00%"* ]] ; then
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
		percents_count=$(echo "${header}" | tr -cd '%' | wc -c)
		echo "${header}"
		Genus=$(echo "${header}" | cut -d' ' -f1 | cut -d'-' -f3)
		species=$(echo "${header}" | cut -d' ' -f2 | cut -d'(' -f1 | sed 's/[][]//g')
		confidence_index=$(echo "${header}" | cut -d' ' -f1 | cut -d'-' -f1,2)
		echo "${Genus}-${species}"
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

# Check if genus was assigned
if [[ ! -z ${Genus} ]]; then
	Genus=$(echo ${Genus} | tr -d [:space:] | tr -d "[]")
fi
# Check if species was assigned
if [[ ! -z ${species} ]]; then
	species=$(echo ${species} | tr -d [:space:])
fi

# Check if genus was assigned as peptoclostridium and relabel it as Clostridium for downstream analyses relying on this older naming convention
if [[ ${Genus} == "Peptoclostridium" ]]; then
	Genus="Clostridium"
fi


if [[ -f ${tax_DB} ]]; then
	# Using premade database fill in upper levels of taxonomy info based on genus
	while IFS= read -r line  || [ -n "$line" ]; do
		DB_genus=$(echo ${line} | cut -d"," -f1)
		#echo ":${Genus}:${DB_genus}:"
		if [[ "${Genus,}" = "${DB_genus}" ]]; then
				Domain=$(echo "${line}" | cut -d"," -f2)
				Phylum=$(echo "${line}" | cut -d"," -f3)
				Class=$(echo "${line}" | cut -d"," -f4)
				Order=$(echo "${line}" | cut -d"," -f5)
				Family=$(echo "${line}" | cut -d"," -f6 | tr -d '\r' )
				#echo ":${Family}:"
				break
		fi
	done < "${tax_DB}"
else
	echo "taxes.csv (${tax_DB}) does not exist"
fi

# Print output to tax file for sample
#echo -e "(${source})-${confidence_index}-${source_file}\nD:	${Domain}\nP:	${Phylum}\nC:	${Class}\nO:	${Order}\nF:	${Family}\nG:	${Genus}\ns:	${species}\n" > "${OUTDATADIR}/${sample_name}.tax"
echo -e "(${source})-${confidence_index}-${source_file}\nD:	${Domain}\nP:	${Phylum}\nC:	${Class}\nO:	${Order}\nF:	${Family}\nG:	${Genus}\ns:	${species}\n" > "${sample_name}.tax"
