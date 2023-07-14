#!/bin/bash -l
#
# Description: Checks sample output folders for correct fiels and tests the thresholds for passability. Edited for use in P
#
# Usage: ./pipeline_stats_writer.sh -e explicit_path_to_isolate_folder
#
# Output location: results/ID/ID.synopsis
#
# Modules required: None
#
# v1.0 (05/25/2022)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
  echo "Usage: pipeline_stats_writer.sh args(* are required)
    -a raw_read_counts.txt
    -b total_read_counts.txt
    -c removed_adapter_R1.fastq.gz
    -d removed_adapter_R2.fastq.gz
    -e kraken2_trimd_report
    -f kraken2_trimd_summary
    -g krona_trimd.html
    "
}

kraken2_unclass_flag=30
kraken2_contamination_threshold=25
ani_coverage_threshold=80

# Parse command line options
options_found=0
while getopts ":1?a:b:c:d:e:f:g:x:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
    a)
      #echo "Option -a triggered, argument = ${OPTARG}"
      raw_read_counts=${OPTARG};;
    b)
      #echo "Option -b triggered, argument = ${OPTARG}"
      total_read_counts=${OPTARG};;
    c)
      #echo "Option -c triggered, argument = ${OPTARG}"
      removed_adapter_R1=${OPTARG};;
    d)
      #echo "Option -d triggered, argument = ${OPTARG}"
      removed_adapter_R2=${OPTARG};;
    e)
      #echo "Option -e triggered, argument = ${OPTARG}"
      kraken2_trimd_report=${OPTARG};;
    f)
      #echo "Option -f triggered, argument = ${OPTARG}"
      kraken2_trimd_summary=${OPTARG};;
    g)
      #echo "Option -g triggered, argument = ${OPTARG}"
      krona_trimd=${OPTARG};;
    x)
      #echo "Option -w triggered, argument = ${OPTARG}"
      srst2_file=${OPTARG};;
    :)
      echo "Option -${OPTARG} requires as argument";;
    1)
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

 # Checks for proper argumentation
 # Checks for correct parameter s and sets appropriate outdatadirs

#OUTDATADIR="${epath}"
#project=$(echo "${epath}" | rev | cut -d'/' -f2 | rev)
sample_name=$(basename "${raw_read_counts}" _raw_read_counts.txt)


# Creates and prints header info for the sample being processed
today=$(date)
echo "----------Checking ${sample_name} for successful completion on ----------"  > "${sample_name}.synopsis"
#echo "Sample output folder starts at: " "${OUTDATADIR}"
status="SUCCESS"

printf "%-30s: %-8s : %s\\n" "Summarized" "SUCCESS" "${today}"  >> "${sample_name}.synopsis"

# Crude check to see if any early steps were performed on reads to indicate this was a full run or assembly only
if [[ "${assembly_only}" = "true" ]]; then
  run_type="assembly-only"
else
  run_type="all"
fi


if [[ "${run_type}" == "all" ]]; then
  #Checking existence of FASTQ files
  raw_length_R1=-3
  raw_length_R2=-3
  #Checking QC counts
  if [[ -s "${raw_read_counts}" ]]; then
    raw_length_R1=$(tail -n1  "${raw_read_counts}" | cut -d$'\t' -f3)
    raw_length_R2=$(tail -n1  "${raw_read_counts}" | cut -d$'\t' -f5)
    if [[ "${raw_length_R1}" -gt 0 ]] && [[ "${raw_length_R2}" -gt 0 ]]; then
      printf "%-30s: %-8s : %s\\n" "FASTQs" "SUCCESS" "R1: ${raw_length_R1}bps R2: ${raw_length_R2}bps"  >> "${sample_name}.synopsis"
    else
      if [[ "${raw_length_R1}" -le 0 ]]; then
        printf "%-30s: %-8s : %s\\n" "FASTQs_R1" "FAILED" "R1 is not represented correctly in reads QC file (${raw_read_counts})"  >> "${sample_name}.synopsis"
        status="FAILED"
        raw_length_R1=0
      else
        printf "%-30s: %-8s : %s\\n" "FASTQs_R1" "SUCCESS" "${raw_length_R1}bps"  >> "${sample_name}.synopsis"
      fi
      if [[ "${raw_length_R2}" -le 0 ]]; then
        printf "%-30s: %-8s : %s\\n" "FASTQs_R2" "FAILED" "R2 is not represented correctly in reads QC file (${raw_read_counts})"  >> "${sample_name}.synopsis"
        status="FAILED"
        raw_length_R2=0
      else
        printf "%-30s: %-8s : %s\\n" "FASTQs_R2" "SUCCESS" "${raw_length_R2}bps"  >> "${sample_name}.synopsis"
      fi
    fi
    raw_exists="true"
  else
    printf "%-30s: %-8s : %s\\n" "FASTQs" "FAILED" "${raw_read_counts} reads QC file does not exist"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "RAW_READ_COUNTS" "FAILED" "${raw_read_counts} reads QC file does not exist"  >> "${sample_name}.synopsis"
    status="FAILED"
    raw_length_R1=0
    raw_length_R2=0
    raw_exists="false"
  fi

  if [[ "${raw_exists}" = "true" ]]; then
    raw_reads=$(tail -n1  "${raw_read_counts}" | cut -d$'\t' -f17)
    raw_pairs=$((raw_reads/2))
    raw_Q30_R1=$(tail -n1 "${raw_read_counts}" | cut -d$'\t' -f14)
    raw_Q30_R1_rounded=$(echo "${raw_Q30_R1}"  | cut -d'.' -f2)
    raw_Q30_R1_rounded=$(echo "${raw_Q30_R1_rounded::2}")
    raw_Q30_R2=$(tail -n1 "${raw_read_counts}" | cut -d$'\t' -f15)
    raw_Q30_R2_rounded=$(echo "${raw_Q30_R2}"  | cut -d'.' -f2)
    raw_Q30_R2_rounded=$(echo "${raw_Q30_R2_rounded::2}")
  if [[ "${raw_reads}" -le 1000000 ]] && [[ "${raw_reads}" -ge 1 ]]; then
      printf "%-30s: %-8s : %s\\n" "RAW_READ_COUNTS" "WARNING" "Low individual read count before trimming: ${raw_reads} (${raw_pairs} paired reads)"  >> "${sample_name}.synopsis"
      status="WARNING"
    elif [[ "${raw_reads}" -le 0 ]]; then
      printf "%-30s: %-8s : %s\\n" "RAW_READ_COUNTS" "FAILED" "No individual read count before trimming: ${raw_reads} (${raw_pairs} paired reads)"  >> "${sample_name}.synopsis"
      status="FAILED"
    else
      printf "%-30s: %-8s : %s\\n" "RAW_READ_COUNTS" "SUCCESS" "${raw_reads} individual reads found in sample (${raw_pairs} paired reads)"  >> "${sample_name}.synopsis"
    fi
    if [[ "${raw_Q30_R1_rounded}" -lt 90 ]]; then
      printf "%-30s: %-8s : %s\\n" "RAW_Q30_R1%" "WARNING" "Q30_R1% at ${raw_Q30_R1_rounded}% (Threshold is 90%)"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
        status="WARNING"
      fi
    else
      printf "%-30s: %-8s : %s\\n" "RAW_Q30_R1%" "SUCCESS" "Q30_R1% at ${raw_Q30_R1_rounded}% (Threshold is 90)"  >> "${sample_name}.synopsis"
    fi
    if [[ "${raw_Q30_R2_rounded}" -lt 70 ]]; then
      printf "%-30s: %-8s : %s\\n" "RAW_Q30_R2%" "WARNING" "Q30_R2% at ${raw_Q30_R2_rounded}% (Threshold is 70)"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
        status="WARNING"
      fi
    else
      printf "%-30s: %-8s : %s\\n" "RAW_Q30_R2%" "SUCCESS" "Q30_R2% at ${raw_Q30_R2_rounded}% (Threshold is 70)"  >> "${sample_name}.synopsis"
    fi
    trimmed_exists="true"
  else
    printf "%-30s: %-8s : %s\\n" "RAW_READ_COUNTS" "FAILED" "${sample_name}_raw_read_counts.txt not found"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "RAW_Q30_R1%" "FAILED" "${sample_name}_raw_read_counts.txt not found"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "RAW_Q30_R2%" "FAILED" "${sample_name}_raw_read_counts.txt not found"  >> "${sample_name}.synopsis"
    status="FAILED"
    trimmed_exists="false"
  fi


  #Checking existence of FASTQ files
  trimmed_length_R1=-3
  trimmed_length_R2=-3
  trimmed_length_unpaired=-3
  #Checking QC counts
  if [[ -s "${total_read_counts}" ]]; then
    trimmed_length_R1=$(tail -n1  ${total_read_counts} | cut -d$'\t' -f3)
    trimmed_length_R2=$(tail -n1  ${total_read_counts} | cut -d$'\t' -f5)
    trimmed_length_unpaired=$(tail -n1  ${total_read_counts} | cut -d$'\t' -f7)
    if [[ "${trimmed_length_R1}" -gt 0 ]] && [[ "${trimmed_length_R2}" -gt 0 ]] && [[ "${trimmed_length_unpaired}" -gt 0 ]]; then
      printf "%-30s: %-8s : %s\\n" "TRIMMED_BPS" "SUCCESS" "R1: ${trimmed_length_R1}bps R2: ${trimmed_length_R2}bps Unpaired: ${trimmed_length_unpaired}bps"  >> "${sample_name}.synopsis"
      #printf "%-30s: %-8s : %s\\n" "TRIMMED_Unpaired" "SUCCESS" "${trimmed_length_unpaired}bps"
    else
      if [[ "${trimmed_length_R1}" -le 0 ]]; then
        printf "%-30s: %-8s : %s\\n" "TRIMMED_R1" "FAILED" "R1 is not represented correctly in reads QC file (${total_read_counts})"  >> "${sample_name}.synopsis"
        status="FAILED"
        trimmed_length_R1=0
      else
        printf "%-30s: %-8s : %s\\n" "TRIMMED_R1" "SUCCESS" "${trimmed_length_R1}bps"  >> "${sample_name}.synopsis"
      fi
      if [[ "${trimmed_length_R2}" -le 0 ]]; then
        printf "%-30s: %-8s : %s\\n" "TRIMMED_R2" "FAILED" "R2 is not represented correctly in reads QC file (${total_read_counts})"  >> "${sample_name}.synopsis"
        status="FAILED"
        trimmed_length_R2=0
      else
        printf "%-30s: %-8s : %s\\n" "TRIMMED_R2" "SUCCESS" "${trimmed_length_R2}bps"  >> "${sample_name}.synopsis"
      fi
      if [[ "${trimmed_length_unpaired}" -gt 0 ]]; then
        printf "%-30s: %-8s : %s\\n" "TRIMMED_UNPAIRED" "SUCCESS" "${trimmed_length_unpaired}bps"  >> "${sample_name}.synopsis"
      else
        printf "%-30s: %-8s : %s\\n" "TRIMMED_UNPAIRED" "ALERT" "No orphaned reads were found"  >> "${sample_name}.synopsis"
        trimmed_length_unpaired=0
      fi
    fi
    total_exists="true"
  else
    printf "%-30s: %-8s : %s\\n" "TRIMMED_FASTQs" "FAILED" "${trimmed_read_counts} reads QC file does not exist"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "TRIMMED_READ_COUNTS" "FAILED" "${trimmed_read_counts} reads QC file does not exist"  >> "${sample_name}.synopsis"
    status="FAILED"
    trimmed_length_R1=0
    trimmed_length_R2=0
    trimmed_length_unpaired=0
    total_exists="false"
  fi

  total_trimmed_reads=-3

  if [[ "${total_exists}" ]]; then
    total_trimmed_reads=$(tail -n1  "${total_read_counts}" | cut -d$'\t' -f23)
    orphaned_reads=$(tail -n1  "${total_read_counts}" | cut -d$'\t' -f6)
    trimmed_reads=$(( total_trimmed_reads - orphaned_reads))
    paired_trimmed=$((trimmed_reads/2))

    trimmed_Q30_R1=$(tail -n1 "${total_read_counts}" | cut -d'	' -f19)
    trimmed_Q30_R1_rounded=$(echo "${trimmed_Q30_R1}"  | cut -d'.' -f2)
    trimmed_Q30_R1_rounded=$(echo "${trimmed_Q30_R1_rounded::2}")
    trimmed_Q30_R2=$(tail -n1 "${total_read_counts}" | cut -d'	' -f20)
    trimmed_Q30_R2_rounded=$(echo "${trimmed_Q30_R2}"  | cut -d'.' -f2)
    trimmed_Q30_R2_rounded=$(echo "${trimmed_Q30_R2_rounded::2}")
    trimmed_Q30_unpaired=$(tail -n1 "${total_read_counts}" | cut -d'	' -f21)
    trimmed_Q30_unpaired_rounded=$(echo "${trimmed_Q30_unpaired}"  | cut -d'.' -f2)
    trimmed_Q30_unpaired_rounded=$(echo "${trimmed_Q30_unpaired_rounded::2}")
    if [[ "${total_trimmed_reads}" -le 1000000 ]] && [[ "${total_trimmed_reads}" -ge 1 ]]; then
      printf "%-30s: %-8s : %s\\n" "TRIMMED_READ_COUNTS" "WARNING" "Low individual read count after trimming: ${total_trimmed_reads} (${paired_trimmed} paired reads, ${orphaned_reads} singled reads)"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "ALERT" ]] || [[ "${status}" = "SUCCESS" ]]; then
        status="WARNING"
      fi
    elif [[ "${total_trimmed_reads}" -le 0 ]]; then
      printf "%-30s: %-8s : %s\\n" "TRIMMED_READ_COUNTS" "FAILED" "No individual read count after trimming: ${trimmed_reads} (${paired_trimmed} paired reads, ${orphaned_reads} singled reads)"  >> "${sample_name}.synopsis"
      status="FAILED"
    else
      printf "%-30s: %-8s : %s\\n" "TRIMMED_READ_COUNTS" "SUCCESS" "${total_trimmed_reads} individual reads found in sample (${paired_trimmed} paired reads, ${orphaned_reads} singled reads)"  >> "${sample_name}.synopsis"
    fi
    if [[ "${trimmed_Q30_R1_rounded}" -lt 90 ]]; then
      printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R1%" "WARNING" "Q30_R1% at ${trimmed_Q30_R1_rounded}% (Threshold is 90%)"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
        status="WARNING"
      fi
    else
      printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R1%" "SUCCESS" "Q30_R1% at ${trimmed_Q30_R1_rounded}% (Threshold is 90)"  >> "${sample_name}.synopsis"
    fi
    if [[ "${trimmed_Q30_R2_rounded}" -lt 70 ]]; then
      printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R2%" "WARNING" "Q30_R2% at ${trimmed_Q30_R2_rounded}% (Threshold is 70)"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
        status="WARNING"
      fi
    else
      printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R2%" "SUCCESS" "Q30_R2% at ${trimmed_Q30_R2_rounded}% (Threshold is 70)"  >> "${sample_name}.synopsis"
    fi
  else
    printf "%-30s: %-8s : %s\\n" "TRIMMED_READ_COUNTS" "FAILED" "${sample_name}_trimmed_read_counts.txt not found"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R1%" "FAILED" "${sample_name}_trimmed_read_counts.txt not found"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R2%" "FAILED" "${sample_name}_trimmed_read_counts.txt not found"  >> "${sample_name}.synopsis"
    status="FAILED"
  fi


	# #Checking fastP output folder
	# remAdapt_length_R1=-1
	# if [[ -s "${removed_adapter_R1}" ]]; then
  #   remAdapt_length_R1=$(zcat ${removed_adapter_R1} | paste - - - - | cut -f2 | tr -d '\n' | wc -c)
	# 	remAdapt_R1_diff=$(( raw_length_R1 - remAdapt_length_R1 ))
	# 	if [[ ${raw_length_R1} -gt 0 ]]; then
	# 		if [[ "${remAdapt_length_R1}" -gt 0 ]]; then
 	# 			R1_adapt_percent_loss=$(( remAdapt_R1_diff * 100 / ${raw_length_R1} ))
	# 			printf "%-30s: %-8s : %s\\n" "BBDUK-R1" "SUCCESS" "R1: ${remAdapt_length_R1}bps (${R1_adapt_percent_loss}% loss)"  >> "${sample_name}.synopsis"
	# 		else
	# 			printf "%-30s: %-8s : %s\\n" "BBDUK-R1" "WARNING" "R1 FASTQ has no remaining bases"  >> "${sample_name}.synopsis"
	# 			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "SUCCESS" ]]; then
	# 				status="WARNING"
	# 			fi
	# 		fi
	# 	else
	# 		if [[ "${remAdapt_length_R1}" -gt 0 ]]; then
	# 			printf "%-30s: %-8s : %s\\n" "BBDUK-R1" "ALERT" "No raw FASTQs availble. R1: ${remAdapt_length_R1}bps (UNK% loss)"  >> "${sample_name}.synopsis"
	# 			if [[ "${status}" = "SUCCESS" ]]; then
	# 				status="ALERT"
	# 			fi
	# 		else
	# 			printf "%-30s: %-8s : %s\\n" "BBDUK-R1" "WARNING" "No raw FASTQs available. BBDUK FASTQ R1 has no bases"  >> "${sample_name}.synopsis"
	# 			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "SUCCESS" ]]; then
	# 				status="WARNING"
	# 			fi
	# 		fi
	# 	fi
	# else
	# 	printf "%-30s: %-8s : %s\\n" "BBDUK-R1" "FAILED" "No BBDUK R1 FASTQ file found"  >> "${sample_name}.synopsis"
	# fi
  #
  # remAdapt_length_R2=-1
  # if [[ -s "${removed_adapter_R2}" ]]; then
	# 	remAdapt_length_R2=$(zcat ${removed_adapter_R2} | paste - - - - | cut -f2 | tr -d '\n' | wc -c)
	# 	remAdapt_R2_diff=$(( raw_length_R2 - remAdapt_length_R2 ))
	# 	if [[ ${raw_length_R2} -gt 0 ]]; then
	# 		if [[ "${remAdapt_length_R2}" -gt 0 ]]; then
 	# 			R2_adapt_percent_loss=$(( remAdapt_R2_diff * 100 / ${raw_length_R2} ))
	# 			printf "%-30s: %-8s : %s\\n" "BBDUK-R2" "SUCCESS" "R2: ${remAdapt_length_R2}bps (${R2_adapt_percent_loss}% loss)"  >> "${sample_name}.synopsis"
	# 		else
	# 			printf "%-30s: %-8s : %s\\n" "BBDUK-R2" "WARNING" "R2 FASTQ has no remaining bases"  >> "${sample_name}.synopsis"
	# 			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "SUCCESS" ]]; then
	# 				status="WARNING"
	# 			fi
	# 		fi
	# 	else
	# 		if [[ "${remAdapt_length_R2}" -gt 0 ]]; then
	# 			printf "%-30s: %-8s : %s\\n" "BBDUK-R2" "ALERT" "No raw FASTQs availble. R2: ${remAdapt_length_R1}bps (UNK% loss)"  >> "${sample_name}.synopsis"
	# 			if [[ "${status}" = "SUCCESS" ]]; then
	# 				status="ALERT"
	# 			fi
	# 		else
	# 			printf "%-30s: %-8s : %s\\n" "BBDUK-R1" "WARNING" "No raw FASTQs available. BBDUK FASTQ R2 has no bases"  >> "${sample_name}.synopsis"
	# 			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "SUCCESS" ]]; then
	# 				status="WARNING"
	# 			fi
	# 		fi
	# 	fi
	# else
	# 	printf "%-30s: %-8s : %s\\n" "BBDUK-R2" "FAILED" "No BBDUK R2 FASTQ file found"  >> "${sample_name}.synopsis"
	# fi



	#Check kraken2 on preAssembly
	kraken2_pre_success="false"
	if [[ -s "${kraken2_trimd_report}" ]]; then
		#printf "%-30s: %-8s : %s\\n" "kraken2 preassembly" "SUCCESS" "Found"
		kraken2_pre_success="true"
	else
		printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS" "FAILED" "${sample_name}.kraken2_trimd.report.txt not found"  >> "${sample_name}.synopsis"
		status="FAILED"
	fi

	#Check Krona output
	if [[ "${kraken2_pre_success}" = "true" ]]; then
		if [[ -s "${krona_trimd}" ]]; then
			#printf "%-30s: %-8s : %s\\n" "krona-kraken2-preasm" "SUCCESS" "Found"
			:
		else
			printf "%-30s: %-8s : %s\\n" "KRONA_READS" "FAILED" "${sample_name}_trimd.html not found"  >> "${sample_name}.synopsis"
			status="FAILED"
		fi
	else
		printf "%-30s: %-8s : %s\\n" "KRONA_READS" "FAILED" "kraken2 reads did not complete successfully"  >> "${sample_name}.synopsis"
		status="FAILED"
	fi

	#Check extraction and unclassified value
	if [[ -s "${kraken2_trimd_summary}" ]]; then
		# Extracts many elements of the summary file to report unclassified and species classified reads and percentages
		unclass=$(head -n 1 "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)
		#true_unclass=$(head -n 1 "${OUTDATADIR}/kraken/preAssembly/${sample_name}_kraken_summary_paired.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
		domain=$(sed -n '2p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)
		genuspre=$(sed -n '7p' "${kraken2_trimd_summary}" | cut -d' ' -f3 | xargs echo)
		speciespre=$(sed -n '8p' "${kraken2_trimd_summary}" | cut -d' ' -f3- | xargs echo)
		speciesprepercent=$(sed -n '8p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)
    if [[ "${unclass}" = "UNK" ]]; then
      unclass=0
      unclass_string="UNKNOWN"
    else
      unclass_string="${unclass}"
    fi
		#true_speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken/preAssembly/${sample_name}_kraken_summary_paired.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
		# If there are no reads at the domain level, then report no classified reads

    if (( $(echo $domain 0 | awk '{if ($1 <= $2) print 1;}') )); then
		#if (( $(echo "${domain} <= 0" | bc -l) )); then
      if [[ "${kraken2_pre_success}" = true ]]; then
        printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "FAILED" "There are no classified reads"  >> "${sample_name}.synopsis"
        status="FAILED"
      else
        printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "FAILED" "KRAKEN2_READS did not complete successfully"  >> "${sample_name}.synopsis"
        status="FAILED"
      fi
		# If there are classified reads then check to see if percent unclassifed falls above the threshold limit. Report warning if too high or success and stats if below
		else
      if (( $(echo ${unclass} ${kraken2_unclass_flag} | awk '{if ($1 > $2) print 1;}') )); then
			#if (( $(echo "${unclass} > ${kraken2_unclass_flag}" | bc -l) )); then
				printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "WARNING" "unclassified reads comprise ${unclass_string}% of total"  >> "${sample_name}.synopsis"
				if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
					status="WARNING"
				fi
			else
				#printf "%-30s: %-8s : %s\\n" "kraken on reads" "SUCCESS" "${speciespercent}%${true_speciespercent%} ${genuspre} ${speciespre} with ${unclass}%${true_unclass%} unclassified reads"
				printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "SUCCESS" "${speciesprepercent}% ${genuspre} ${speciespre} with ${unclass_string}% unclassified reads"  >> "${sample_name}.synopsis"
			fi
		fi
	# If no summary file was found
	else
		printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "FAILED" "${sample_name}.kraken2_trimd.classifiedreads.txt not found"  >> "${sample_name}.synopsis"
		status="FAILED"
	fi

	# Quick separate check for contamination by finding # of species above ${contamination_threshold} in list file from kraken2
	if [[ -s "${kraken2_trimd_report}" ]]; then
		number_of_species=0
		while IFS= read -r line; do
			arrLine=(${line})
			# First element in array is the percent of reads identified as the current taxa
			percent=${arrLine[0]}
			percent_integer=$(echo "${percent}" | cut -d'.' -f1)
			# 3rd element is the taxon level classification
			# echo "${percent_integer} vs ${contamination_threshold}"
			classification=${arrLine[3]}
			if [[ "${classification}" == "S" ]] && (( percent_integer > kraken2_contamination_threshold )); then
				#echo "Adding ${arrLine[5]}-${percent_integer}-${contamination_threshold} to list"
				number_of_species=$(( number_of_species + 1 ))
			fi
		done < "${kraken2_trimd_report}"

		if [[ "${number_of_species}" -gt 1 ]]; then
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS_CONTAM" "WARNING" "${number_of_species} species have been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
			if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
					status="WARNING"
			fi
		elif [[ "${number_of_species}" -eq 1 ]]; then
			:
		else
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS_CONTAM" "WARNING" "No species have been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
        status="WARNING"
      fi
		fi
		#echo "Number of species: ${number_of_species}"
	fi

else
	printf "%-30s: %-8s : %s\\n" "QC_COUNTS" "NA" "Assembly only isolate"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "Q30_STATS" "NA" "Assembly only isolate"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "BBDUK" "NA" "Assembly only isolate"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "TRIMMING" "NA" "Assembly only isolate"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS" "NA" "Assembly only isolate"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRONA_READS" "NA" "Assembly only isolate"  >> "${sample_name}.synopsis"
fi

#Check spades assembly
if [[ -s "${SPAdes_assembly}" ]]; then
  # Count the number of '>' in the assembly file before trimming
  full_scaffolds=$(zgrep -c '>' ${SPAdes_assembly})
  printf "%-30s: %-8s : %s\\n" "ASSEMBLY" "SUCCESS" "${full_scaffolds} scaffolds found"  >> "${sample_name}.synopsis"
else
  printf "%-30s: %-8s : %s\\n" "ASSEMBLY" "FAILED" "${sample_name}.scaffolds.fa.gz not found"  >> "${sample_name}.synopsis"
  QC_FAIL=$QC_FAIL"smaller_than_1000000_bps(0)-"
  status="FAILED"
fi

kraken2_contamination_threshold=25
printf "%-30s: %-8s : %s\\n" "SRST2" "UNKNOWN" "SPAdes failed ${sample_name}__fullgenes__*.txt was not looked for."  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "SCAFFOLD_TRIM" "FAILED" "${sample_name}.filtered.scaffolds.fa.gz not found"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD" "FAILED" "${sample_name}.kraken2_asmbld.report.txt not found"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRONA_ASMBLD" "FAILED" "kraken2 unweighted did not complete successfully - SPAdes failure"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "FAILED" "kraken2 assembly did not comlete successfully - SPAdes failure"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD_CONTAM" "ALERT" "No Genera have been found above ${kraken2_contamination_threshold}% abundance"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED" "FAILED" "${sample_name}.kraken2_wtasmbld.report.txt not found"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRONA_WEIGHTED" "FAILED" "kraken2 weighted did not complete successfully - SPAdes failure"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "FAILED" "${sample_name}.kraken2_wtasmbld.classifiedreads.txt not found"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD_CONTAM" "ALERT" "No Genera have been found above ${kraken2_contamination_threshold}% abundance"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "QUAST" "FAILED" "${sample_name}_report.tsv does not exist"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "TAXA-${tax_source}" "FAILED" "No Taxa File found"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "ASSEMBLY_RATIO(SD)" "FAILED" "No Ratio File exists"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "COVERAGE" "FAILED" "${avg_coverage}x coverage based on trimmed reads (Min:30x)"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "BUSCO" "FAILED" "short_summary.*.${sample_name}.scaffolds.fa.txt not found"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "FAILED" "NO ${sample_name}_REFSEQ_*.fastANI.txt  file"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "MLST" "FAILED" "${sample_name}.tsv does not exist"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "GAMMA_AR" "FAILED" "${sample_name}_ResGANNCBI_*.gamma does not exist"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "PLASMID_REPLICONS" "FAILED" "${sample_name}_PF-*.gamma does not exist"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "HYPERVIRULENCE" "FAILED" "${sample_name}_Hyper_Virulence_*.gamma does not exist"  >> "${sample_name}.synopsis"
status="FAILED"
QC_FAIL=$QC_FAIL"coverage_below_30($avg_coverage), SPAdes failed"

if [[ -n "${QC_FAIL}" ]]; then
  QC_FAIL=${QC_FAIL%?}
  printf "%-30s: %-8s : %s\\n" "Auto Pass/FAIL" "FAIL" "$QC_FAIL"  >> "${sample_name}.synopsis"
  status="FAILED"
else
  printf "%-30s: %-8s : %s\\n" "Auto Pass/FAIL" "PASS" "Minimum Requirements met for coverage(30x)/ratio_stdev(<2.58)/min_length(>1000000) to pass auto QC filtering"  >> "${sample_name}.synopsis"
fi

echo "---------- ${sample_name} completed as ${status} ----------"  >> "${sample_name}.synopsis"
echo "WARNINGS: out of line with what is expected and MAY cause problems downstream."  >> "${sample_name}.synopsis"
echo "ALERT: something to note, does not mean it is a poor-quality assembly."  >> "${sample_name}.synopsis"


#Script exited gracefully (unless something else inside failed)
exit 0
