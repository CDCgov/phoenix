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
# Created by Nick Vlachos (nvx4@cdc.gov)

#version=2.0 # (11/15/2023) Changed to signify adoption of CLIA minded versioning. This version is equivalent to previous version 1.0.3 (04/19/2022)
version=2.1 # (04/10/2024) update to handle scientific notation

#  Function to print out help blurb
show_help () {
	echo "Usage is ./sort_and_prep_dists.sh  -a assembly -x dists_file -d database_of_fastas_matching_dist_file_output [-V show version]"
	echo "Output is saved to folder where .list file exists"
}

# Parse command line options
options_found=0
while getopts ":h?x:d:a:o:t:V" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
			show_help
			exit 0
			;;
		x)
			echo "Option -x triggered, argument = ${OPTARG}"
			dist_file=${OPTARG}
			;;
		a)
			echo "Option -a triggered, argument = ${OPTARG}"
			assembly_file=${OPTARG}
			;;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			outdir=${OPTARG}
			;;
		t)
			echo "Option -t triggered"
			terra=${OPTARG}
			;;
		V)
			show_version="True"
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
	echo "No argument supplied to sort_and_prep.sh, exiting"
	show_help
	exit 1
fi

if [[ "${show_version}" = "True" ]]; then
	echo "sort_and_prep_dist.sh: ${version}"
	exit
fi

# set the correct path for bc/wget - needed for terra
if [[ $terra = "terra" ]]; then
	wget_path=/opt/conda/envs/phoenix/bin/wget
	python_path=/opt/conda/envs/phoenix/bin/python3
	certificate_check="--no-check-certificate"
else
	wget_path=/usr/bin/wget
	python_path=python3
	certificate_check=""
fi

# Based upon standard naming protocols pulling last portion of path off should result in proper name
sample_name=$(basename "${dist_file}" .txt)

sorted_dists="${dist_file//.txt/_sorted.txt}"

#swapped -n for -g to handle scientific notation
sort -k3 -g -o "${sorted_dists}" "${dist_file}"

cutoff=$(head -n20 "${sorted_dists}" | tail -n1 | cut -d'	' -f3)

### Could add some logic here to prevent terrible isolats from continuing on

echo "Cutoff IS: ${cutoff}"

matches=0
# Needed a new variable to put a hard stop on fill-ins being 150% of orignal. Example - if max ani samples to use is 20, the 30th sample is the last one that could be used as filler
# Temporarily setting value here, should move to main parameters in the future?
max_hits=40 #Target is 20, but we'll allow twice as many even though this likely only occurs in crummy isolates

#echo "${assembly_file}" > "${sample_name}_best_MASH_hits.txt"

## create empty log file to add to 
touch ${sample_name}_genome_download_log.txt

while IFS= read -r var; do
	echo "${var}"
	source=$(echo "${var}" | cut -d$'\t' -f1)
	dist=$(echo ${var} | cut -d' ' -f3)
	kmers=$(echo ${var} | cut -d' ' -f5 | cut -d'/' -f1)
	echo "dist-${dist} - ${source}"
	#use python to handle cases of scientific notation
	result=$($python_path -c "print(1 if float('$dist') <= float('$cutoff') else 0)")
	# Also setting a minimum kmer threshold to ensure 1000 crappy hits dont make it ner the top with 1/1000 kmer matches
	if ([ "$result" -eq 1 ] && [ ${kmers} -gt 5 ] && [ ${matches} -le ${max_hits} ]); then
		if [[ -f "${outdir}/${source}.gz" ]]; then
			echo "${outdir}/${source}.gz" >> "${sample_name}_best_MASH_hits.txt"
			matches=$(( matches + 1))
		else
			# Loop through field numbers 3-6 to find valid GCF_ pattern
			GCF_found=0
			for field_num in {3..6}; do
				GCF_name=$(echo "${source}" | cut -d'_' -f${field_num}-)
				GCF_check=${GCF_name:0:4}
				if [[ "${GCF_check}" = "GCF_" ]]; then
					GCF_found=1
					filename=$(echo ${source} | cut -d'_' -f${field_num}- | rev | cut -d'_' -f2,3,4 | rev)
					alpha=${filename:4:3}
					beta=${filename:7:3}
					charlie=${filename:10:3}
					echo "Copying - ${filename}"
					echo "Trying - wget $certificate_check https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${GCF_name}.gz -O ${outdir}/${source}.gz"
					$wget_path $certificate_check https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${GCF_name}.gz -O ${outdir}/${source}.gz
		#			echo "Trying - wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz -O ${filename}_genomic.fna.gz"
		#			wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz -O ${filename}_genomic.fna.gz
					#curl --remote-name --remote-time "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz"
					# Check if file exists and is not empty, if so it will move on to the next best hit. Dont really know a better way to deal with a crappy download situation
					if [[ -s "${outdir}/${source}.gz" ]]; then
						echo "${outdir}/${source}.gz" >> "${sample_name}_best_MASH_hits.txt"
			#			echo "${GCF_name}.gz" >> "${sample_name}_best_MASH_hits.txt"
						matches=$(( matches + 1))
					else
						echo "${source} did not download correctly" >> ${sample_name}_genome_download_log.txt
					fi
					break  # Exit the for loop once GCF_ is found
				fi
			done
			if [[ ${GCF_found} -eq 0 ]]; then
				echo "GCF check did not pass for any field position in ${source}" >> ${sample_name}_genome_download_log.txt
			fi
		fi
	else
		break
	fi
	counter=$(( counter + 1 ))
done < ${sorted_dists}
