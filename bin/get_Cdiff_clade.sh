#!/bin/bash -l

#$ -o run_GAMA.out
#$ -e run_GAMA.err
#$ -N run_GAMA
#$ -cwd
#$ -q short.q

#
# Description: A very simple tool to pull out the clade of a Cdiff sample from the pubmlst profile lists
#
# Usage: ./run_GAMA_for_normal.sh -e explicit_path_to_isolate_folder [-c config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/MLST/
#
# Modules required: None
#
# v1.1 (09/27/2023)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#


#  Function to print out help blurb
show_help () {
	echo "Usage: ./get_Cdiff_clade.sh -m mlst_combined_file -r path_to_mlst_database_folder"
}

version="1.1"

# Parse command line options
options_found=0
while getopts ":h?m:r:V" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
			show_help
			exit 0;;
		m)
			echo "Option -m triggered, argument = ${OPTARG}"
			input=${OPTARG};;
		r)
			echo "Option -r triggered, argument = ${OPTARG}"
			mlst_db=${OPTARG};;
		V)
			show_version="True";;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

if [[ "${show_version}" = "True" ]]; then
	echo "${version}"
	exit
fi

# Checks for correct parameter s and sets appropriate outdatadirs
if [[ ! -f "${input}" ]]; then
echo "Empty sample path supplied (-e), exiting"
	exit
fi
sample_name=$(echo "${input}" | rev | cut -d'/' -f1 | cut -d'_' -f2- | rev)

# Need sorter for possible PHX outcome styles
# Process for a single line
line_count=$(cat "${input}" | wc -l)
if [[ "${line_count}" -eq 2 ]]; then
	line=$(tail -n1 "${input}")
	type=$(echo "${line}" | cut -d$'\t' -f5)
	if [[ "${type}" = "A-SUB" ]] || [[ "${type}" = "P-SUB" ]] || [[ "${type}" = "novel_allele" ]] || [[ "${type}" = "novel_profile" ]]; then
		type="MLST_type_not_defined"
	fi
# If there is only a header or nothing at at all
elif [[ "${line_Count}" -lt 2 ]]; then
	type="MLST_type_not_found"
# Sort through multiple possible lines to select the best and give look up THAT clade
else
	type="MLST_type_not_found"
	type_source="none"
	while IFS= read -r line; do
		type_source_temp=$(echo "${line}" | cut -d$'\t' -f2)
		type_temp=$(echo "${line}" | cut -d$'\t' -f5)
		if [[ "${type_source_temp}" = "assembly"* ]]; then
			# Check if the current type has a source other than reads. If none set as current and any assembly reset as current since it 
			if [[ "${type_source}" = "none" ]] || [[ "${type_source}" = "assembly" ]]; then
				if [[ "${type_temp}" = "A-SUB" ]] || [[ "${type_temp}" = "P-SUB" ]] || [[ "${type_temp}" = "novel_allele" ]] || [[ "${type_temp}" = "novel_profile" ]]; then
					type="MLST_type_not_defined"
					type_source="assembly"
				else
					type="${type_temp}"
					break
				fi
			# Else see if the reads answer is better than the assembly answer
			elif [[ "${type_source}" = "reads" ]]; then
				if [[ "${type_temp}" = "A-SUB" ]] || [[ "${type_temp}" = "P-SUB" ]] || [[ "${type_temp}" = "novel_allele" ]] || [[ "${type_temp}" = "novel_profile" ]]; then
					if [[ "${type}" = "MLST_type_not_defined" ]]; then
						type="MLST_type_not_defined"
						type_source="assembly"
					# Do nothing if the reads type is a real answer when this assembly one is not
					else
						:
					fi
				else
					type="${type_temp}"
					break
				fi
			# Means it is an assembly source, and any assembly source that hasnt exited is an undefined..do nothing
			else
				:
			fi
		# Do it all again for reads
		elif [[ "${type_source_temp}" = "reads" ]]; then
			# Check if the current type has a source other than reads. If none set as current and any assembly reset as current since it 
			if [[ "${type_source}" = "none" ]]; then
				if [[ "${type_temp}" = "A-SUB" ]] || [[ "${type_temp}" = "P-SUB" ]] || [[ "${type_temp}" = "novel_allele" ]] || [[ "${type_temp}" = "novel_profile" ]]; then
					type="MLST_type_not_defined"
					type_source="reads"
				else
					type="${type_temp}"
					type_source="reads"
				fi
			elif [[ "${type_source}" = "reads" ]]; then
				if [[ "${type_temp}" = "A-SUB" ]] || [[ "${type_temp}" = "P-SUB" ]] || [[ "${type_temp}" = "novel_allele" ]] || [[ "${type_temp}" = "novel_profile" ]]; then
					if [[ "${type}" = "A-SUB" ]] || [[ "${type}" = "P-SUB" ]] || [[ "${type}" = "novel_allele" ]] || [[ "${type}" = "novel_profile" ]]; then
						type="MLST_type_not_defined"
						type_source="reads"
					# Do nothing if there is already a value set based on reads
					else
						:
					fi
				else
					if [[ "${type}" = "A-SUB" ]] || [[ "${type}" = "P-SUB" ]] || [[ "${type}" = "novel_allele" ]] || [[ "${type}" = "novel_profile" ]]; then
						type="${type_temp}"
						type_source="reads"
					# Do nothing if there is already a value set based on reads
					else
						:
					fi
				fi
			# Else see if the reads answer is better than the assembly answer
			elif [[ "${type_source}" = "assembly" ]]; then
				if [[ "${type_temp}" = "A-SUB" ]] || [[ "${type_temp}" = "P-SUB" ]] || [[ "${type_temp}" = "novel_allele" ]] || [[ "${type_temp}" = "novel_profile" ]]; then
					if [[ "${type}" = "MLST_type_not_defined" ]]; then
						type="MLST_type_not_defined"
						type_source="assembly"
					# Do nothing if the reads type is a real answer when this assembly one is not
					else
						:
					fi
				else
					if [[ "${type}" = "A-SUB" ]] || [[ "${type}" = "P-SUB" ]] || [[ "${type}" = "novel_allele" ]] || [[ "${type}" = "novel_profile" ]]; then
						type="${type_temp}"
						type_source="reads"
					# Do nothing if there is already a value set based on reads
					else
						:
					fi
				fi
			fi
		# Dont know what else it could be but adding for catch all
		else	
			:
		fi
	done < "${input}"
fi

echo "Type: ${type}"
profile_line=$(grep -P "^${type}\t" "${mlst_db}"/pubmlst/cdifficile/cdifficile.txt)
echo "Line: ${profile_line}"
if [[ -n "%{profile_line}" ]] && [[ "${profile_line}" != "" ]]; then
    clade=$(echo "${profile_line}" | cut -d$'\t' -f9)
fi
if [[ -z "${clade}" ]] || [[ "${clade}" = '' ]]; then
    clade="clade_Not_defined"
fi

echo -e "Sample name\tClade" > "${sample_name}_cdifficile_clade.tsv"
echo -e "${sample_name}\t${clade}" >> "${sample_name}_cdifficile_clade.tsv"
#echo -e "Sample name\tClade"
#echo -e "${sample_name}\t${clade}"

#Script exited gracefully (unless something else inside failed)
exit 0
