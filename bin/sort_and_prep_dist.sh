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
	echo "Usage is ./sort_and_prep_dists.sh  -a assembly -x dists_file -d database_of_fastas_matching_dist_file_output"
	echo "Output is saved to folder where .list file exists"
}

# Parse command line options
options_found=0
while getopts ":h?x:d:a:t:" option; do
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
	echo "No argument supplied to best_hit_from_kraken_noconfig.sh, exiting"
	show_help
	exit 1
fi

# set the correct path for bc/wget - needed for terra
if [[ $terra = "terra" ]]; then
	bc_path=/opt/conda/envs/phoenix/bin/bc
	wget_path=/opt/conda/envs/phoenix/bin/wget
	certificate_check="--no-check-certificate"
else
	bc_path=bc
	wget_path=/usr/bin/wget
	certificate_check=""
fi

# Based upon standard naming protocols pulling last portion of path off should result in proper name
sample_name=$(basename "${dist_file}" .txt)

sorted_dists="${dist_file//.txt/_sorted.txt}"

sort -k3 -n -o "${sorted_dists}" "${dist_file}"

cutoff=$(head -n20 "${sorted_dists}" | tail -n1 | cut -d'	' -f3)

### Could add some logic here to prevent terrible isolats from continuing on

echo "Cutoff IS: ${cutoff}"

matches=0
# Needed a new variable to put a hard stop on fill-ins being 150% of orignal. Example - if max ani samples to use is 20, the 30th sample is the last one that could be used as filler
#maxplus_ani_samples=$(echo 'scale=0; 1.5 * '${max_ani_samples}' / 1' | bc -l)

#echo "${assembly_file}" > "${sample_name}_best_MASH_hits.txt"

while IFS= read -r var; do
	echo "${var}"
	source=$(echo "${var}" | cut -d$'\t' -f1)
	dist=$(echo ${var} | cut -d' ' -f3)
	kmers=$(echo ${var} | cut -d' ' -f5 | cut -d'/' -f1)
	echo "dist-${dist} - ${source}"
		if ((( $(echo "$dist <= $cutoff" | $bc_path -l) )) && [ ${kmers} -gt 0 ]); then
		if [[ -f "${source}.gz" ]]; then
			echo "${source}.gz" >> "${sample_name}_best_MASH_hits.txt"
#		if [[ -f "${GCF_name}.gz" ]]; then
#			echo "${GCF_name}.gz" >> "${sample_name}_best_MASH_hits.txt"
			matches=$(( matches + 1))
		else
			filename=$(echo ${source} | cut -d'_' -f3- | rev | cut -d'_' -f2,3,4 | rev)
			GCF_name=$(echo "${source}" | cut -d'_' -f3-)
			GCF_check=${GCF_name:0:4}
			if [[ "${GCF_check}" = "GCF_" ]]; then
				alpha=${filename:4:3}
				beta=${filename:7:3}
				charlie=${filename:10:3}
				echo "Copying - ${filename}"
				echo "Trying - wget $certificate_check https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${GCF_name}.gz -O ${source}.gz"
				$wget_path $certificate_check https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${GCF_name}.gz -O ${source}.gz
	#			echo "Trying - wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz -O ${filename}_genomic.fna.gz"
	#			wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz -O ${filename}_genomic.fna.gz
				#curl --remote-name --remote-time "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz"
				echo "${source}.gz" >> "${sample_name}_best_MASH_hits.txt"
	#			echo "${GCF_name}.gz" >> "${sample_name}_best_MASH_hits.txt"
				matches=$(( matches + 1))
			else
				echo "GCF check did not pass, look into the differences of ${source}"
			fi
		fi
	else
		break
	fi
	counter=$(( counter + 1 ))
done < ${sorted_dists}
