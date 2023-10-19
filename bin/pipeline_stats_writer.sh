#!/bin/bash -l
#
# Description: Checks sample output folders for correct fiels and tests the thresholds for passability. Edited for use in PHoeNIx
#
# Usage: ./pipeline_stats_writer.sh -e explicit_path_to_isolate_folder
#
# Output location: results/ID/ID.synopsis
#
# Modules required: None
#
# v1.1 (03/22/2023)
#
# Created by Nick Vlachos (nvx4@cdc.gov), edits by Jill Hagey (qpk9@cdc.gov)
#

#  Function to print out help blurb
show_help () {
  echo "Usage: pipeline_stats_writer.sh args(* are required)
    -a raw_read_counts.txt
    -b total_read_counts.txt
    -c gc_content_file
    -d sample_name
    -e kraken2_trimd_report
    -f kraken2_trimd_summary
    -g krona_trimd.html
    -h* assembly
    -i* trimmed_assembly
    -j* kraken2_asmbld_report
    -k* kraken2_asmbled_summary
    -l* krona_asmbld.html
    -m* kraken2_weighted_report
    -n* kraken2_weighted_summary
    -o* krona_weighted.html
    -p* quast_report.tsv
    -q* taxID.tax
    -r* assembly_ratio_file
    -s* busco_specific_short_summary
    -t* fastANI_formatted_file
    -u* gamma_AR.gamma
    -v* gamma_replicon.gamma
    -w* gamma_HV.gamma
    -x srst_fullgenes_file
    -y* MLST_file (or more)
    -z* assembly_only_sample (true or false)
    -1* amr_tsv_file
    -5 coverage
    -3
    "
}

kraken2_unclass_flag=30
kraken2_contamination_threshold=25
ani_coverage_threshold=80

# Parse command line options
options_found=0
while getopts ":1?a:b:c:d:e:f:g:h:i:j:k:l:m:n:o:p:q:r:s:t:u:v:w:x:y:z:4:2:5:3" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
    a)
      echo "Option -a triggered, argument = ${OPTARG}"
      raw_read_counts=${OPTARG};;
    b)
      echo "Option -b triggered, argument = ${OPTARG}"
      total_read_counts=${OPTARG};;
    c)
      echo "Option -c triggered, argument = ${OPTARG}"
      gc_content_file=${OPTARG};;
    d)
      echo "Option -d triggered, argument = ${OPTARG}"
      sample_name=${OPTARG};;
    e)
      echo "Option -e triggered, argument = ${OPTARG}"
      kraken2_trimd_report=${OPTARG};;
    f)
      echo "Option -f triggered, argument = ${OPTARG}"
      kraken2_trimd_summary=${OPTARG};;
    g)
      echo "Option -g triggered, argument = ${OPTARG}"
      krona_trimd=${OPTARG};;
    h)
      echo "Option -h triggered, argument = ${OPTARG}"
      SPAdes_assembly=${OPTARG};;
    i)
      echo "Option -i triggered, argument = ${OPTARG}"
      trimmed_assembly=${OPTARG};;
    j)
      echo "Option -j triggered, argument = ${OPTARG}"
      kraken2_asmbld_report=${OPTARG};;
    k)
      echo "Option -k triggered, argument = ${OPTARG}"
      kraken2_asmbled_summary=${OPTARG};;
    l)
      echo "Option -l triggered, argument = ${OPTARG}"
      krona_asmbld=${OPTARG};;
    m)
      echo "Option -m triggered, argument = ${OPTARG}"
      kraken2_weighted_report=${OPTARG};;
    n)
      echo "Option -n triggered, argument = ${OPTARG}"
      kraken2_weighted_summary=${OPTARG};;
    o)
      echo "Option -o triggered, argument = ${OPTARG}"
      krona_weighted=${OPTARG};;
    p)
      echo "Option -p triggered, argument = ${OPTARG}"
      quast_report=${OPTARG};;
    q)
      echo "Option -q triggered, argument = ${OPTARG}"
      taxID_file=${OPTARG};;
    r)
      echo "Option -r triggered, argument = ${OPTARG}"
      assembly_ratio_file=${OPTARG};;
    s)
      echo "Option -s triggered, argument = ${OPTARG}"
      busco_summary=${OPTARG};;
    t)
      echo "Option -t triggered, argument = ${OPTARG}"
      formatted_fastANI=${OPTARG};;
    u)
      echo "Option -t triggered, argument = ${OPTARG}"
      gamma_AR=${OPTARG};;
    v)
      echo "Option -u triggered, argument = ${OPTARG}"
      gamma_replicon=${OPTARG};;
    w)
      echo "Option -v triggered, argument = ${OPTARG}"
      gamma_HV=${OPTARG};;
    x)
      echo "Option -w triggered, argument = ${OPTARG}"
      srst2_file=${OPTARG};;
    y)
      echo "Option -x triggered, argument = ${OPTARG}"
      mlst_file=${OPTARG};;
    z)
      echo "Option -z triggered, argument = ${OPTARG}"
      assembly_only="true";;
    4)
      echo "Option -1 triggered, argument = ${OPTARG}"
      amr_file=${OPTARG};;
		2)
			echo "Option -2 triggered"
			terra=${OPTARG};;
    5)
      echo "Option -5 triggered, argument = ${OPTARG}"
      coverage=${OPTARG};;
    3)
      echo "Option -3 triggered, argument = ${OPTARG}"
      internal_phoenix="true";;
    :)
      echo "Option -${OPTARG} requires as argument";;
    1)
      show_help
      exit 0
      ;;
	esac
done

# set the correct path for bc/wget - needed for terra
if [[ $terra = "terra" ]]; then
	bc_path=/opt/conda/envs/phoenix/bin/bc
else
	bc_path=bc
fi

#check coverage parameter if >30 use that number
if [ "${coverage}" -ge 30 ]; then
  reads_min=${coverage}
	reads_low=40
	reads_high=150
else
  reads_min=30
	reads_low=40
	reads_high=150
fi

if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

 # Checks for proper argumentation
 # Checks for correct parameter s and sets appropriate outdatadirs

#OUTDATADIR="${epath}"
#project=$(echo "${epath}" | rev | cut -d'/' -f2 | rev)
#sample_name=$(basename "${raw_read_counts}" _raw_read_counts.txt) # don't need this anymore as we are passing with -1


# Creates and prints header info for the sample being processed
today=$(date)
echo "---------- Checking ${sample_name} for successful completion on ${today} ----------"  > "${sample_name}.synopsis"
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
      printf "%-30s: %-8s : %s\\n" "RAW_Q30_R1%" "SUCCESS" "Q30_R1% at ${raw_Q30_R1_rounded}% (Threshold is 90%)"  >> "${sample_name}.synopsis"
    fi
    if [[ "${raw_Q30_R2_rounded}" -lt 70 ]]; then
      printf "%-30s: %-8s : %s\\n" "RAW_Q30_R2%" "WARNING" "Q30_R2% at ${raw_Q30_R2_rounded}% (Threshold is 70%)"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
        status="WARNING"
      fi
    else
      printf "%-30s: %-8s : %s\\n" "RAW_Q30_R2%" "SUCCESS" "Q30_R2% at ${raw_Q30_R2_rounded}% (Threshold is 70%)"  >> "${sample_name}.synopsis"
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

  if [[ "${total_exists}" = true ]]; then
    total_trimmed_reads=$(tail -n1  "${total_read_counts}" | cut -d$'\t' -f24)
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
      printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R1%" "SUCCESS" "Q30_R1% at ${trimmed_Q30_R1_rounded}% (Threshold is 90%)"  >> "${sample_name}.synopsis"
    fi
    if [[ "${trimmed_Q30_R2_rounded}" -lt 70 ]]; then
      printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R2%" "WARNING" "Q30_R2% at ${trimmed_Q30_R2_rounded}% (Threshold is 70%)"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
        status="WARNING"
      fi
    else
      printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R2%" "SUCCESS" "Q30_R2% at ${trimmed_Q30_R2_rounded}% (Threshold is 70%)"  >> "${sample_name}.synopsis"
    fi
  else
    printf "%-30s: %-8s : %s\\n" "TRIMMED_READ_COUNTS" "FAILED" "${sample_name}_trimmed_read_counts.txt not found"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R1%" "FAILED" "${sample_name}_trimmed_read_counts.txt not found"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R2%" "FAILED" "${sample_name}_trimmed_read_counts.txt not found"  >> "${sample_name}.synopsis"
    status="FAILED"
  fi

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
		unclass=$(sed -n '2p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)
		domain=$(sed -n '3p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)
		genuspre=$(sed -n '8p' "${kraken2_trimd_summary}" | cut -d' ' -f3 | xargs echo)
		speciespre=$(sed -n '9p' "${kraken2_trimd_summary}" | cut -d' ' -f3- | xargs echo)
		speciesprepercent=$(sed -n '9p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)
        genusprepercent=$(sed -n '8p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)

    ### Prepping for header addition
    #unclass=$(sed -n '2p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)
		#domain=$(sed -n '3p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)
		#genuspre=$(sed -n '8p' "${kraken2_trimd_summary}" | cut -d' ' -f3 | xargs echo)
		#speciespre=$(sed -n '9p' "${kraken2_trimd_summary}" | cut -d' ' -f3- | xargs echo)
		#speciesprepercent=$(sed -n '10p' "${kraken2_trimd_summary}" | cut -d' ' -f2 | xargs echo)

    if [[ "${unclass}" = "UNK" ]]; then
      unclass=0
      unclass_string="UNKNOWN"
    else
      unclass_string="${unclass}"
    fi
		#true_speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken/preAssembly/${sample_name}_kraken_summary_paired.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
		# If there are no reads at the domain level, then report no classified reads

    #if (( $(echo $domain 0 | awk '{if ($1 <= $2) print 1;}') )); then
		if (( $(echo "${domain} <= 0" | $bc_path -l) )); then
      if [[ "${kraken2_pre_success}" = true ]]; then
        printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "FAILED" "There are no classified reads"  >> "${sample_name}.synopsis"
        status="FAILED"
      else
        printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "FAILED" "KRAKEN2_READS did not complete successfully"  >> "${sample_name}.synopsis"
        status="FAILED"
      fi
		# If there are classified reads then check to see if percent unclassifed falls above the threshold limit. Report warning if too high or success and stats if below
		else
      #if (( $(echo ${unclass} ${kraken2_unclass_flag} | awk '{if ($1 > $2) print 1;}') )); then
			if (( $(echo "${unclass} > ${kraken2_unclass_flag}" | $bc_path -l) )); then
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
		number_of_genera=0
		while IFS= read -r line; do
			arrLine=(${line})
			# First element in array is the percent of reads identified as the current taxa
			percent=${arrLine[0]}
			percent_integer=$(echo "${percent}" | cut -d'.' -f1)
			# 3rd element is the taxon level classification
			# echo "${percent_integer} vs ${contamination_threshold}"
			classification=${arrLine[3]}
			if [[ "${classification}" == "G" ]] && (( percent_integer > kraken2_contamination_threshold )); then
				#echo "Adding ${arrLine[5]}-${percent_integer}-${contamination_threshold} to list"
				number_of_genera=$(( number_of_genera + 1 ))
			fi
		done < "${kraken2_trimd_report}"

		if [[ "${number_of_genera}" -gt 1 ]]; then
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS_CONTAM" "WARNING" "${number_of_genera} genera have been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
			if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
					status="WARNING"
			fi
		elif [[ "${number_of_genera}" -eq 1 ]]; then
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS_CONTAM" "SUCCESS" "Only one genus has been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
		else
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS_CONTAM" "WARNING" "No genera have been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
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



#Check short scaffolds reduction script
if [[ -s "${trimmed_assembly}" ]]; then
	# Count the number of '>' still remaining after trimming the contig file
	full_longies=$(zgrep -c '>' "${trimmed_assembly}")
	# Calculate the number of lost (short) scaffolds
	full_shorties=$(( full_scaffolds - full_longies ))
	if [ -z ${full_shorties} ]; then
		full_shorties=0
	fi
	#echo "${full_longies}"
	if [[ "${full_longies}" -le 200 ]]; then
		printf "%-30s: %-8s : %s\\n" "SCAFFOLD_TRIM" "SUCCESS" "${full_longies} scaffolds remain. ${full_shorties} were removed due to shortness"  >> "${sample_name}.synopsis"
	else
		printf "%-30s: %-8s : %s\\n" "SCAFFOLD_TRIM" "WARNING" "${full_longies} scaffolds remain which is high. ${full_shorties} were removed due to shortness"  >> "${sample_name}.synopsis"
		if [[ "${status}" == "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
			status="WARNING"
		fi
	fi
else
	printf "%-30s: %-8s : %s\\n" "SCAFFOLD_TRIM" "FAILED" "${sample_name}.filtered.scaffolds.fa.gz not found"  >> "${sample_name}.synopsis"
	status="FAILED"
fi

if [[ "${internal_phoenix}" == "true" ]]; then
	#Check kraken2 on assembly
	kraken2_unweighted_success="false"
	if [[ -s "${kraken2_asmbld_report}" ]]; then
		#printf "%-30s: %-8s : %s\\n" "kraken2 postassembly" "SUCCESS" "Found"
		kraken2_unweighted_success="true"
	else
		printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD" "FAILED" "${sample_name}.kraken2_asmbld.report.txt not found"  >> "${sample_name}.synopsis"
		status="FAILED"
	fi

	#Check Krona output of assembly
	if [[ "${kraken2_unweighted_success}" = "true" ]]; then
		if [[ -s "${krona_asmbld}" ]]; then
			#printf "%-30s: %-8s : %s\\n" "krona-kraken2-pstasm" "SUCCESS" "Found"
			:
		else
			printf "%-30s: %-8s : %s\\n" "KRONA_ASMBLD" "FAILED" "${sample_name}_asmbld.html not found"  >> "${sample_name}.synopsis"
			status="FAILED"
		fi
	else
		printf "%-30s: %-8s : %s\\n" "KRONA_ASMBLD" "FAILED" "kraken2 unweighted did not complete successfully"  >> "${sample_name}.synopsis"
		status="FAILED"
	fi

	#Check extraction and unclassified values for kraken2 post assembly
	if [[ -s "${kraken2_asmbled_summary}" ]]; then
		# Extracts many elements of the summary file to report unclassified and species classified reads and percentages
		unclass=$(sed -n '2p' "${kraken2_asmbled_summary}" | cut -d' ' -f2 | xargs echo)
		domain=$(sed -n '3p' "${kraken2_asmbled_summary}" | cut -d' ' -f2 | xargs echo)
		genuspost=$(sed -n '8p' "${kraken2_asmbled_summary}" | cut -d' ' -f3 | xargs echo)
		speciespost=$(sed -n '9p' "${kraken2_asmbled_summary}" | cut -d' ' -f3- | xargs echo)
		speciespostpercent=$(sed -n '9p' "${kraken2_asmbled_summary}" | cut -d' ' -f2 | xargs echo)
		genuspostpercent=$(sed -n '8p' "${kraken2_asmbled_summary}" | cut -d' ' -f2 | xargs echo)

    #unclass=$(sed -n '2p' "${kraken2_asmbled_summary}" | cut -d' ' -f2 | xargs echo)
		#domain=$(sed -n '3p' "${kraken2_asmbled_summary}" | cut -d' ' -f2 | xargs echo)
		#genuspost=$(sed -n '8p' "${kraken2_asmbled_summary}" | cut -d' ' -f3 | xargs echo)
		#speciespost=$(sed -n '9p' "${kraken2_asmbled_summary}" | cut -d' ' -f3- | xargs echo)
		#speciespostpercent=$(sed -n '9p' "${kraken2_asmbled_summary}" | cut -d' ' -f2 | xargs echo)
		#genuspostpercent=$(sed -n '8p' "${kraken2_asmbled_summary}" | cut -d' ' -f2 | xargs echo)



		if [[ "${unclass}" = "UNK" ]]; then
			unclass=0
			unclass_string="UNKNOWN"
		else
			unclass_string="${unclass}"
		fi
		#true_speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken2/postAssembly/${sample_name}_kraken_summary_assembled.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
		# If there are no reads at the domain level, then report no classified reads

		#if (( $(echo $domain 0 | awk '{if ($1 <= $2) print 1;}') )); then
		if (( $(echo "${domain} <= 0" | $bc_path -l) )); then
			if [[ "${kraken2_unweighted_success}" = true ]]; then
				printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "FAILED" "There are no classified reads (Did post assembly kraken2 fail too?)"	>> "${sample_name}.synopsis"
				status="FAILED"
			else
				printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "FAILED" "kraken2 assembly did not complete successfully"	>> "${sample_name}.synopsis"
				status="FAILED"
			fi
		# If there are classified reads then check to see if percent unclassifed falls above the threshold limit. Report warning if too high or success and stats if below
		else
			#if (( $(echo ${unclass} ${kraken2_unclass_flag} | awk '{if ($1 > $2) print 1;}') )); then
			if (( $(echo "${unclass} > ${kraken2_unclass_flag}" | $bc_path -l) )); then
				printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "WARNING" "unclassified scaffolds comprise ${unclass_string}% of total"	>> "${sample_name}.synopsis"
				if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
					status="WARNING"
				fi
			#elif (( $(echo ${genuspostpercent} 50 | awk '{if ($1 < $2) print 1;}') )); then
			elif (( $(echo "${genuspostpercent} < 50" | $bc_path -l) )); then
				printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "WARNING" "Genus-${genuspost}(${genuspostpercent}%) is under 50% (species ${speciespost} (${speciespostpercent}%)), possibly contaminated or scaffolds are weighted unevenly"  >> "${sample_name}.synopsis"
				if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
					status="WARNING"
				fi
			else
				#printf "%-30s: %-8s : %s\\n" "kraken2 on assembly" "SUCCESS" "${speciespercent}%${true_speciespercent%} ${genuspost} ${speciespost} with ${unclass}%${true_unclass%} unclassified contigs"
				printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "SUCCESS" "${genuspost}(${genuspostpercent}%) ${speciespost}(${speciespostpercent}%) with ${unclass_string}% unclassified scaffolds"  >> "${sample_name}.synopsis"
			fi
		fi
	# If no summary file was found
	else
		printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "FAILED" "${sample_name}.kraken2_asmbld.classifiedreads.txt not found"  >> "${sample_name}.synopsis"
		status="FAILED"
	fi

	# Quick separate check for contamination by finding # of species above ${contamination_threshold} in list file from kraken2
	if [[ -s "${kraken2_asmbld_report}" ]]; then
		number_of_genera=0
		#echo "${kraken2_asmbld_report}" >> "${sample_name}.synopsis"
		while IFS= read -r line; do
			arrLine=(${line})
			# First element in array is the percent of reads identified as the current taxa
			percent=${arrLine[0]}
  		percent_integer=$(echo "${percent}" | cut -d'.' -f1)
				# 3rd element is the taxon level classification
				classification=${arrLine[3]}
				#printf "${line} - ${classification} - ${percent_integer} - ${kraken2_contamination_threshold}" >> "${sample_name}.synopsis"
				if [[ "${classification}" = "G" ]] && (( percent_integer > kraken2_contamination_threshold )); then
					number_of_genera=$(( number_of_genera + 1 ))
					#printf "adding ${classification} at ${percent_integer} (above ${kraken2_contamination_threshold})" >> "${sample_name}.synopsis"
				fi
		done < "${kraken2_asmbld_report}"
		###echo "${number_of_genera}"
		if [[ $number_of_genera -gt 1 ]]; then
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD_CONTAM" "ALERT" "${number_of_genera} genera have been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
			if [[ "${status}" == "SUCCESS" ]]; then
				status="ALERT"
			fi
		elif [[ ${number_of_genera} -eq 1 ]]; then
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD_CONTAM" "SUCCESS" "Only one genus has been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
		else
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD_CONTAM" "ALERT" "No genera have been found above ${kraken2_contamination_threshold}% abundance"  >> "${sample_name}.synopsis"
			if [[ "${status}" = "ALERT" ]] || [[ "${status}" = "SUCCESS" ]]; then
				status="WARNING"
			fi
			#echo "Number of genera: ${number_of_genera}"
		fi
  fi
fi

kraken2_weighted_success=false
if [[ -s "${kraken2_weighted_report}" ]]; then
	#printf "%-30s: %-8s : %s\\n" "kraken2 weighted" "SUCCESS" "Found"
	kraken2_weighted_success=true
else
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED" "FAILED" "${sample_name}.kraken2_wtasmbld.report.txt not found"  >> "${sample_name}.synopsis"
	status="FAILED"
fi

#Check Krona output of weighted assembly
if [[ "${kraken2_weighted_success}" = true ]]; then
  if [[ -s "${krona_weighted}" ]]; then
    #printf "%-30s: %-8s : %s\\n" "krona-kraken2-pstasm" "SUCCESS" "Found"
    :
  else
    printf "%-30s: %-8s : %s\\n" "KRONA_WEIGHTED" "FAILED" "${sample_name}_weighted.html not found"  >> "${sample_name}.synopsis"
    status="FAILED"
  fi
else
  printf "%-30s: %-8s : %s\\n" "KRONA_WEIGHTED" "FAILED" "kraken2 weighted did not complete successfully"  >> "${sample_name}.synopsis"
  status="FAILED"
fi

#Check extraction and unclassified values for weighted kraken2 post assembly
if [[ -s "${kraken2_weighted_summary}" ]]; then
	# Extracts many elements of the summary file to report unclassified and species classified reads and percentages
	unclass=$(sed -n '2p' "${kraken2_weighted_summary}" | cut -d' ' -f2 | xargs echo)
	domain=$(sed -n '3p' "${kraken2_weighted_summary}" | cut -d' ' -f2 | xargs echo)
	genusweighted=$(sed -n '8p' "${kraken2_weighted_summary}" | cut -d' ' -f3 | xargs echo)
	speciesweighted=$(sed -n '9p' "${kraken2_weighted_summary}" | cut -d' ' -f3- | xargs echo)
	speciesweightedpercent=$(sed -n '9p' "${kraken2_weighted_summary}" | cut -d' ' -f2 | xargs echo)
	genusweightedpercent=$(sed -n '8p' "${kraken2_weighted_summary}" | cut -d' ' -f2 | xargs echo)


  #unclass=$(sed -n '2p' "${kraken2_weighted_summary}" | cut -d' ' -f2 | xargs echo)
	#domain=$(sed -n '3p' "${kraken2_weighted_summary}" | cut -d' ' -f2 | xargs echo)
	#genusweighted=$(sed -n '8p' "${kraken2_weighted_summary}" | cut -d' ' -f3 | xargs echo)
	#speciesweighted=$(sed -n '9p' "${kraken2_weighted_summary}" | cut -d' ' -f3- | xargs echo)
	#speciesweightedpercent=$(sed -n '9p' "${kraken2_weighted_summary}" | cut -d' ' -f2 | xargs echo)
	#genusweightedpercent=$(sed -n '8p' "${kraken2_weighted_summary}" | cut -d' ' -f2 | xargs echo)

	if [[ "${unclass}" = "UNK" ]]; then
		unclass=0
		unclass_string="UNKNOWN"
	else
		unclass_string="${unclass}"
	fi
	#true_speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken2/postAssembly/${sample_name}_kraken_summary_assembled_BP.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
	# If there are no reads at the domain level, then report no classified reads
	#if (( $(echo $domain 0 | awk '{if ($1 <= $2) print 1;}') )); then
	if (( $(echo "${domain} <= 0" | $bc_path -l) )); then
		if [[ "${kraken2_weighted_success}" = true ]]; then
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "FAILED" "There are no classified reads"	>> "${sample_name}.synopsis"
			status="FAILED"
		else
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "FAILED" "Kraken2 weighted did not complete successfully"	>> "${sample_name}.synopsis"
			status="FAILED"
		fi
	# If there are classified reads then check to see if percent unclassifed falls above the threshold limit. Report warning if too high or success and stats if below
	else
		#if (( $(echo ${unclass} ${kraken2_unclass_flag} | awk '{if ($1 > $2) print 1;}') )); then
		if (( $(echo "${unclass} > ${kraken2_unclass_flag}" | $bc_path -l) )); then
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "WARNING" "unclassified reads comprise ${unclass_string}% of total"	>> "${sample_name}.synopsis"
			if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
				status="WARNING"
			fi
		#elif (( $(echo ${genusweightedpercent} 50 | awk '{if ($1 < $2) print 1;}') )); then
		elif (( $(echo "${genusweightedpercent} < 50" | $bc_path -l) )); then
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "FAILED" "Genus-${genusweighted} is under 50% (species-${speciesweighted} ${speciesweightedpercent}%), likely contaminated"  >> "${sample_name}.synopsis"
			status="FAILED"
		else
			#printf "%-30s: %-8s : %s\\n" "weighted Classify" "SUCCESS" "${speciespercent}%${true_speciespercent%} ${genusweighted} ${speciesweighted} with ${unclass}%${true_unclass%} unclassified weighted"
			printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "SUCCESS" "${genusweighted}(${genusweightedpercent}%) ${speciesweighted}(${speciesweightedpercent}%) with ${unclass_string}% unclassified scaffolds"  >> "${sample_name}.synopsis"
		fi
	fi
# If no summary file was found
else
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "FAILED" "${sample_name}.kraken2_wtasmbld.classifiedreads.txt not found"  >> "${sample_name}.synopsis"
	status="FAILED"
fi

# Quick separate check for contamination by finding # of species above ${contamination_threshold} in list file from kraken2
if [[ -s "${kraken2_weighted_report}" ]]; then
	number_of_genera=0
	root=0
	unclass=0
	total=0
	while IFS= read -r line; do
		arrLine=(${line})
		# First element in array is the percent of reads identified as the current taxa
		if [[ "${arrLine[3]}" = "U" ]]; then
			unclass=${arrLine[0]}
		elif [[ "${arrLine[3]}" = "R" ]]; then
			root=${arrLine[0]}
			total_percent=$(echo "${unclass} + ${root}" | $bc_path)
			#echo "total percent:${unclass} + ${root} = ${total}"
    else
      percent=$(echo "${arrLine[0]} ${total_percent}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
      #percent=${arrLine[0]}
      percent_integer=$(echo "${percent}" | cut -d'.' -f1)
      # 3rd element is the taxon level classification
      classification=${arrLine[3]}
      #echo "${arrLine[0]} - ${total_percent} - ${percent} - ${percent_integer} - ${kraken2_contamination_threshold} - ${classification}"
      if [[ "${classification}" == "G" ]] && (( percent_integer > kraken2_contamination_threshold )); then
        number_of_genera=$(( number_of_genera + 1 ))
        #echo "adding ${classification} at ${percent_integer} (above ${kraken2_contamination_threshold})"
      fi
    fi
	done < "${kraken2_weighted_report}"

	if [[ $number_of_genera -gt 1 ]]; then
		printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED_CONTAM" "WARNING" "${number_of_genera} genera have been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
		if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
      status="WARNING"
    fi
	elif [[ "${number_of_genera}" -eq 1 ]]; then
		printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED_CONTAM" "SUCCESS" "Only one genus has been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
	else
		printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED_CONTAM" "WARNING" "No genera have been found above the ${kraken2_contamination_threshold}% threshold"  >> "${sample_name}.synopsis"
    if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
      status="WARNING"
    fi
	fi
	#echo "Number of species: ${number_of_genera}"
fi

#Check QUAST
if [[ -s "${quast_report}" ]]; then
	# Extract the useful bits and report (to compare to Toms)
	contig_num=$(sed -n '14p' "${quast_report}" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
	assembly_length=$(sed -n '16p' "${quast_report}" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
	N50=$(sed -n '18p' "${quast_report}" | sed -r 's/[\t]+/ /g'| cut -d' ' -f2 )
	GC_con=$(sed -n '17p' "${quast_report}" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
	#contig_num=$(sed -n '14p' "${quast_report}" | rev | cut -d$'\t ' -f1 | rev)
	#assembly_length=$(sed -n '16p' "${quast_report}" | rev | cut -d$'\t ' -f1 | rev)
	#N50=$(sed -n '18p' "${quast_report}" | rev | cut -d$'\t ' -f1 | rev)
	#GC_con=$(sed -n '17p' "${quast_report}" | rev | cut -d$'\t ' -f1 | rev)
	printf "%-30s: %-8s : %s\\n" "QUAST" "SUCCESS" "#-${contig_num} length-${assembly_length} n50-${N50} %GC-${GC_con}"  >> "${sample_name}.synopsis"
	if [[ "${assembly_length}" -lt 1000000 ]]; then
		QC_FAIL=$QC_FAIL"smaller_than_1000000_bps(${assembly_length})-"
	fi
else
	printf "%-30s: %-8s : %s\\n" "QUAST" "FAILED" "${sample_name}_report.tsv does not exist"  >> "${sample_name}.synopsis"
	status="FAILED"
fi

# Getting GC% information out of quast and gc files
if [[ -s "${quast_report}" ]]; then
	# Extract gc content from quast report
	GC_con=$(sed -n '17p' "${quast_report}" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3 )
  if [[ -f "${gc_content_file}" ]]; then # if gc_content file exists get stdev line
    gc_stdev=$(awk 'NR==3' "${gc_content_file}" | cut -d' ' -f2)
    #checks if "Not calculated on species with n<10 references"
    if [ "${gc_stdev}" = "Not" ]; then #check that there was enough references for stdev calculation
      printf "%-30s: %-8s : %s\\n" "QUAST_GC_Content" "ALERT" "Low References for STDev - ${GC_con}x(${gc_stdev}-SD)"  >> "${sample_name}.synopsis"
    else
      #Check that GC content is within the range we want (2.58 on either side)
      gc_mean=$(awk 'NR==6' "${gc_content_file}" | cut -d' ' -f2)
      devs=$( echo 2.58 \* ${gc_stdev} | $bc_path -l)
      gc_right=$( echo ${gc_mean} + ${devs} | $bc_path -l) # right of the mean
      gc_left=$( echo ${gc_mean} - ${devs} | $bc_path -l) # left of the mean
      # if "Not calculated on species with n<10 references" in stdev 
      if [[ $gc_stdev == "Not" ]]; then
        gc_stdev=0
      else
        gc_stdev=$( printf "%.5f" ${gc_stdev})
      fi
      if [[ $GC_con < $gc_left ]]; then
        gc_left=$( echo ${gc_mean} - ${devs} | $bc_path -l | xargs printf "%.5f") # just making it more readable
        gc_stdev=$( printf "%.5f" ${gc_stdev}) # just making it more readable
        gc_mean=$( printf "%.5f" ${gc_mean}) # just making it more readable
        printf "%-30s: %-8s : %s\\n" "QUAST_GC_Content" "WARNING" "GC Content (%GC-${GC_con}) is less than ${gc_left} (2.58*${gc_stdev}stdevs) away from the mean of ${gc_mean}."  >> "${sample_name}.synopsis"
      elif [[ $GC_con > $gc_right ]]; then
        gc_right=$( echo ${gc_mean} + ${devs} | $bc_path -l | xargs printf "%.5f") # just making it more readable
        gc_stdev=$( printf "%.5f" ${gc_stdev}) # just making it more readable
        gc_mean=$( printf "%.5f" ${gc_mean}) # just making it more readable
        printf "%-30s: %-8s : %s\\n" "QUAST_GC_Content" "WARNING" "GC Content (%GC-${GC_con}) is greater than ${gc_right} (2.58*${gc_stdev}stdevs) away from the mean of ${gc_mean}."  >> "${sample_name}.synopsis"
      else
        gc_right=$( echo ${gc_mean} + ${devs} | $bc_path -l | xargs printf "%.5f") # just making it more readable
        gc_left=$( echo ${gc_mean} - ${devs} | $bc_path -l | xargs printf "%.5f") # just making it more readable
        gc_mean=$( printf "%.5f" ${gc_mean}) # just making it more readable
        printf "%-30s: %-8s : %s\\n" "QUAST_GC_Content" "SUCCESS" "%GC-${GC_con} is within ${gc_left}-${gc_right} (2.58*${gc_stdev}stdevs) away from the mean of ${gc_mean}."  >> "${sample_name}.synopsis"
      fi
    fi
  else # if gc_content file does not exist report failure
    printf "%-30s: %-8s : %s\\n" "QUAST_GC_Content" "FAILED" "%GC-${GC_con}, but ${sample_name}_GC_content_20210819.txt does not exists. ${sample_name}_GC_content_20210819.txt was suppose to be created in assembly ratio step."  >> "${sample_name}.synopsis"
  fi
else
  printf "%-30s: %-8s : %s\\n" "QUAST_GC_Content" "FAILED" "${sample_name}_report.tsv does not exist"  >> "${sample_name}.synopsis"
	status="FAILED"
fi

if [[ -s "${taxID_file}" ]]; then
	source_call=$(head -n1 "${taxID_file}")
	tax_source="UNK"
  tax_source=$(echo "${source_call}" | cut -d$'\t' -f1)
  dec_family=$(grep "F:" "${taxID_file}" | cut -d$'\t' -f2)
  dec_genus=$(grep "G:" "${taxID_file}" | cut -d$'\t' -f2)
  dec_species=$(grep "s:" "${taxID_file}" | cut -d$'\t' -f2)

	if [[ "$dec_genus" != "Not_assigned" ]] && [[ "$dec_species" != "Not_assigned" ]]; then
		printf "%-30s: %-8s : %s\\n" "TAXA-${tax_source}" "SUCCESS" "${dec_genus} ${dec_species}"  >> "${sample_name}.synopsis"
	elif [[ "$dec_genus" != "Not_assigned" ]]; then
		printf "%-30s: %-8s : %s\\n" "TAXA" "FAILED" "None of the classifiers completed successfully"  >> "${sample_name}.synopsis"
	elif [[ "$dec_genus" != "Not_assigned" ]]; then
		printf "%-30s: %-8s : %s\\n" "TAXA" "FAILED" "None of the classifiers completed successfully"  >> "${sample_name}.synopsis"
	elif [[ "$dec_species" != "Not_assigned" ]]; then
		printf "%-30s: %-8s : %s\\n" "TAXA" "WARNING" "No Species was able to be determined"  >> "${sample_name}.synopsis"
	fi
else
	printf "%-30s: %-8s : %s\\n" "TAXA-${tax_source}" "FAILED" "No Taxa File found"  >> "${sample_name}.synopsis"
fi

genus_initial="${dec_genus:0:1}"
assembly_ID="${genus_initial}.${dec_species}"
if [[ -f "${assembly_ratio_file}" ]]; then
  ratio_db_date=$(echo "${assembly_ratio_file}" | rev | cut -d'_' -f1 | rev | cut -d'.' -f1)
  assembly_ratio=$(tail -n1 "${assembly_ratio_file}" | cut -d' ' -f2)
  stdev_line=$(head -n4 "${assembly_ratio_file}" | tail -n1)
  species_stdev_line=$(head -n3 "${assembly_ratio_file}" | tail -n1)
  gc_content=$(head  )
    if [[ "${stdev_line}" = "Isolate_St.Devs: N/A" ]]; then
      st_dev="N/A"
    else
      st_dev=$(head -n4 "${assembly_ratio_file}" | tail -n1 | cut -d' ' -f2)
    fi

    #if (( $(echo $assembly_ratio 0 | awk '{if ($1 < $2) print 1;}') )); then
    if (( $(echo "$assembly_ratio < 0" | $bc_path -l) )); then
      printf "%-30s: %-8s : %s\\n" "ASSEMBLY_RATIO(SD)" "WARNING" "No Reference - ${assembly_ratio}x(${st_dev}-SD) against ${assembly_ID}"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
        status="WARNING"
      fi
    #QC_FAIL=$QC_FAIL"STDev_NOREF-"
  elif [[ "${species_stdev_line}" = *"Not calculated on species with n<10 references"* ]] || [[ "${st_dv}" = "N/A" ]]; then
      printf "%-30s: %-8s : %s\\n" "ASSEMBLY_RATIO(SD)" "ALERT" "Low References for STDev - ${assembly_ratio}x(${st_dev}-SD) against ${assembly_ID}"  >> "${sample_name}.synopsis"
      if [[ "${status}" = "SUCCESS" ]]; then
        status="ALERT"
      fi
    #elif (( $(echo $st_dev 2.58 | awk '{if ($1 > $2) print 1;}') )); then
    elif (( $(echo "$st_dev > 2.58" | $bc_path -l) )); then
      printf "%-30s: %-8s : %s\\n" "ASSEMBLY_RATIO(SD)" "FAILED" "St. dev. too large - ${assembly_ratio}x(${st_dev}-SD) against ${assembly_ID}"  >> "${sample_name}.synopsis"
      status="FAILED"
      QC_FAIL=$QC_FAIL"STDev_above_2.58($st_dev)-"
      #elif (( $(echo "$assembly_ratio < 0.8" | $bc_path -l) )); then
      #	printf "%-30s: %-8s : %s\\n" "Assembly ratio" "FAILED" "Too small - ${assembly_ratio}x(${st_dev}-SD) against ${assembly_ID} (DB up to date! Most current DB: ${NCBI_ratio_date})"
      #	status="FAILED"
    else
      printf "%-30s: %-8s : %s\\n" "ASSEMBLY_RATIO(SD)" "SUCCESS" "${assembly_ratio}x(${st_dev}-SD) against ${assembly_ID}"  >> "${sample_name}.synopsis"
    fi
  else
    printf "%-30s: %-8s : %s\\n" "ASSEMBLY_RATIO(SD)" "FAILED" "No Ratio File exists"  >> "${sample_name}.synopsis"
    if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]] || [[ "${status}" = "WARNING" ]]; then
      status="FAILED"
    fi
  fi

	if [[ -s "${total_read_counts}" ]]; then
		bps_post_all=$(tail -n1 "${total_read_counts}" | cut -d'	' -f22)
	# IFS='	' read -r -a qcs <<< "${line}"
	# read_qc_info=${qcs[@]:1}
	# # Extract q30 reads from qcCounts to calculate average coverage as q30_reads/assembly_length
	# q30_reads=$(echo "${read_qc_info}" | awk -F ' ' '{print $2}')
	# # Change later to AWK as this wont work on ASPEN, but consolidate won't likely be run on cluster
	if [[ ${assembly_length} -gt 0 ]] && [[ ${bps_post_all} -gt 0 ]]; then
		avg_coverage=$($bc_path <<<"scale=2 ; ${bps_post_all} / ${assembly_length}")
	else
		avg_coverage=0
	fi
  #echo "trimmed-${avg_coverage}"
  #if (( $(echo $domain 0 | awk '{if ($1 > $2) print 1;}') )) && (( $(echo ${avg_coverage} ${reads_high} | awk '{if ($1 < $2) print 1;}') )); then
	if (( $(echo "${avg_coverage} > ${reads_low}" | $bc_path -l) )) && (( $(echo "${avg_coverage} < ${reads_high}" | $bc_path -l) )); then
		printf "%-30s: %-8s : %s\\n" "COVERAGE" "SUCCESS" "${avg_coverage}x coverage based on trimmed reads (Target:40x, Cutoff: ${reads_min}x)"  >> "${sample_name}.synopsis"
	#elif (( $(echo ${avg_coverage} ${reads_high} | awk '{if ($1 > $2) print 1;}') )); then
  elif (( $(echo "${avg_coverage} > ${reads_high}" | $bc_path -l) )); then
		printf "%-30s: %-8s : %s\\n" "COVERAGE" "ALERT" "${avg_coverage}x coverage based on trimmed reads (Target:<150x)"  >> "${sample_name}.synopsis"
		if [[ "${status}" == "SUCCESS" ]]; then
			status="ALERT"
		fi
  #elif (( $(echo ${avg_coverage} ${reads_min} | awk '{if ($1 > $2) print 1;}') )); then
  elif (( $(echo "${avg_coverage} > ${reads_min}" | $bc_path -l) )); then
		printf "%-30s: %-8s : %s\\n" "COVERAGE" "ALERT" "${avg_coverage}x coverage based on trimmed reads (Target:40x, Cutoff: ${reads_min}x)"  >> "${sample_name}.synopsis"
		if [[ "${status}" == "SUCCESS" ]]; then
			status="ALERT"
		fi
  #elif (( $(echo ${avg_coverage} ${reads_min} | awk '{if ($1 < $2) print 1;}') )); then
	elif (( $(echo "${avg_coverage} < ${reads_min}" | $bc_path -l) )); then
		printf "%-30s: %-8s : %s\\n" "COVERAGE" "FAILED" "${avg_coverage}x coverage based on trimmed reads (Min:30x)"  >> "${sample_name}.synopsis"
		status="FAILED"
    QC_FAIL=$QC_FAIL"coverage_below_30($avg_coverage)-"
	fi
fi

if [[ "${internal_phoenix}" == "true" ]]; then
  #Check BUSCO
  if [[ -s "${busco_summary}" ]]; then
	  # Reads each line of the busco output file to extract the 3 that contain summary data to report
	  while IFS= read -r line; do
      # If the line contains info for found buscos, total buscos, or database info grab it
      if [[ "${line}" == *"Complete BUSCOs (C)"* ]]; then
        #echo "C-"${line}
			  found_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
  	  elif [[ "${line}" == *"Total BUSCO groups searched"* ]]; then
        #echo "T-"${line}
        total_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
  	  elif [[ "${line}" == *"The lineage dataset is:"* ]]; then
        #echo "L-"${line}
        db=$(echo "${line}" | awk -F ' ' '{print $6}')
      fi
      done < "${busco_summary}"
      BUSCO_organism=$(sed 's/_odb[0-9]\+//' <<< $db)
      percent_BUSCO_present=$($bc_path<<<"${found_buscos}*100/${total_buscos}")
      if [[ "${percent_BUSCO_present}" -gt 90 ]]; then
        printf "%-30s: %-8s : %s\\n" "BUSCO_${db^^}" "SUCCESS" "${percent_BUSCO_present}% of expected core genes for ${BUSCO_organism} found (${found_buscos}/${total_buscos})"  >> "${sample_name}.synopsis"
      else
        printf "%-30s: %-8s : %s\\n" "BUSCO_${db^^}" "FAILED" "only ${percent_BUSCO_present}% of expected core genes for ${BUSCO_organism} found (${found_buscos}/${total_buscos})"  >> "${sample_name}.synopsis"
        status="FAILED"
      fi
    # If the busco summary file does not exist
  else
    printf "%-30s: %-8s : %s\\n" "BUSCO" "FAILED" "short_summary.*.${sample_name}.scaffolds.fa.txt not found"  >> "${sample_name}.synopsis"
    status="FAILED"
  fi
fi

#Check fastANI REFSEQ
if [[ -s "${formatted_fastANI}" ]]; then
  #fastANI_date=$(echo "${formatted_fastANI}" | rev | cut -d'_' -f1 | rev | cut -d'.' -f2)
  #fastANI_info=$(head -n2 "${formatted_fastANI}" | tail -n1)
  #percent_match=$(echo "${fastANI_info}" | cut -d'.' -f1)
  #coverage_match=$(echo "${fastANI_info}" | cut -d'-' -f2 | cut -d'.' -f1)

  ### In anticipation of cleaned up output
  fastANI_date=$(echo "${formatted_fastANI}" | rev | cut -d'_' -f1 | rev | cut -d'.' -f2)
  fastANI_info=$(head -n2 "${formatted_fastANI}" | tail -n1)
  percent_match=$(echo "${fastANI_info}" | cut -d'.' -f1)
  #coverage_match=$(echo "${fastANI_info}" | cut -d'-' -f2 | cut -d'.' -f1)
  coverage_match=$(echo "${fastANI_info}" | cut -d$'\t' -f2 | cut -d'.' -f1)
  organism=$(echo "${fastANI_info}" | cut -d$'\t' -f3)
  reference=$(echo "${fastANI_info}" | cut -d$'\t' -f4)

  if [[ "${percent_match}" = "0."* ]]; then
    printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "FAILED" "No assembly file to work with"  >> "${sample_name}.synopsis"
  else
    if [[ "${percent_match}" -ge 95 ]] && [[ "${coverage_match}" -ge ${ani_coverage_threshold} ]]; then
      printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "SUCCESS" "${percent_match}%ID ${coverage_match}%cov  tax=${organism}  ref=${reference}"  >> "${sample_name}.synopsis"
    else
      if [[ "${percent_match}" -lt 95 ]]; then
        if [[ "${coverage_match}" -lt ${ani_coverage_threshold} ]]; then
          printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "FAILED" "% Identity(${percent_match}%) and % coverage(${coverage_match}%) is too low. ${fastANI_info}"  >> "${sample_name}.synopsis"
        else
          printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "FAILED" "% Identity(${percent_match}%) is too low: ${fastANI_info}"  >> "${sample_name}.synopsis"
        fi
      elif [[ "${coverage_match}" -lt ${ani_coverage_threshold} ]]; then
        printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "FAILED" "% coverage is too low (${coverage_match}%). ${fastANI_info}"  >> "${sample_name}.synopsis"
      fi
    fi
  fi
else
  printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "FAILED" "NO ${sample_name}_REFSEQ_*.fastANI.txt  file"  >> "${sample_name}.synopsis"
  status="FAILED"
fi

# check MLST minimus
if [[ -s "${mlst_file}" ]]; then
  line_count=$(wc -l "${mlst_file}" | cut -d' ' -f1)
  echo "${line_count}"
  current_index=1
  if [[ ${line_count} -ge 2 ]]; then
      while read mlst_line; do
          if [[ "${current_index}" -ge 2 ]]; then
              echo "${mlst_line}"
              mlst_db=$(echo "${mlst_line}" | cut -d$'\t' -f4)
              mlst_type=$(echo "${mlst_line}" | cut -d$'\t' -f5)
              mlst_source=$(echo "${mlst_line}" | cut -d$'\t' -f2)
              if [[ "${mlst_db}" = '-' ]]; then
                printf "%-30s: %-8s : %s\\n" "MLST" "FAILED" "No scheme identified"  >> "${sample_name}.synopsis"
              elif [[ "${mlst_type}" = '-' ]]; then
                printf "%-30s: %-8s : %s\\n" "MLST-${mlst_db^^}" "FAILED" "No type identified, but scheme is ${mlst_db} via ${mlst_source}"  >> "${sample_name}.synopsis"
              else
                printf "%-30s: %-8s : %s\\n" "MLST-${mlst_db^^}" "SUCCESS" "ST${mlst_type} via ${mlst_source}"  >> "${sample_name}.synopsis"
              fi
          fi
          current_index=$(( current_index + 1 ))
      done < "${mlst_file}"
  fi
else
  main_line=$(head -n1 ${mlst_file})
  if [[ "${main_line}" = "Sample  Source"* ]]; then
    printf "%-30s: %-8s : %s\\n" "MLST" "FAILED" "No MLST entries in ${sample_name}.tsv"  >> "${sample_name}.synopsis"
  else
    printf "%-30s: %-8s : %s\\n" "MLST" "FAILED" "${sample_name}.tsv does not exist"  >> "${sample_name}.synopsis"
  fi
  status="FAILED"
fi

#Check GAMMA
gamma_genes="NA"
if [[ -s "${gamma_AR}" ]]; then
  gamma_genes=0
  ARDB=$(echo ${gamma_AR} | rev | cut -d'_' -f2,3 | rev)
  while read line_in; do
    if [[ "${line_in}" = "Gene	Contig	Start	Stop	Match_Type"* ]]; then
      :
    else
      gamma_genes=$(( gamma_genes + 1 ))
    fi
  done < "${gamma_AR}"
  if [[ ${gamma_genes} -eq 0 ]]; then
    printf "%-30s: %-8s : %s\\n" "GAMMA_AR" "SUCCESS" "No AR genes were found from ${ARDB}"  >> "${sample_name}.synopsis"
  elif [[ ${gamma_genes} -ge 1 ]]; then
    printf "%-30s: %-8s : %s\\n" "GAMMA_AR" "SUCCESS" "${gamma_genes} AR gene(s) found from ${ARDB}"  >> "${sample_name}.synopsis"
  else
    echo "Should never get here (gamma_genes less than 0)"
  fi
else
  printf "%-30s: %-8s : %s\\n" "GAMMA_AR" "FAILED" "${sample_name}_ResGANNCBI_*.gamma does not exist"  >> "${sample_name}.synopsis"
  status="FAILED"
fi

#Check NCBI AMRFinder Genes
amr_genes="NA"
echo "file name"
echo ${amr_file}
if [[ -s "${amr_file}" ]]; then
  amr_genes=0
  echo "Got here"
  while read line_in; do
    if [[ "${line_in}" = *"POINT"* ]]; then
      amr_genes=$(( amr_genes + 1 ))
    else
      :
    fi
  done < "${amr_file}"
  if [[ "${amr_genes}" -eq 0 ]]; then
    printf "%-30s: %-8s : %s\\n" "AMRFINDER" "SUCCESS" "No point mutations were found"  >> "${sample_name}.synopsis"
  elif [[ "${amr_genes}" -ge 1 ]]; then
    printf "%-30s: %-8s : %s\\n" "AMRFINDER" "SUCCESS" "${amr_genes} point mutation(s) found"  >> "${sample_name}.synopsis"
  else
    echo "Should never get here (amr_genes less than 0)"
  fi
else
  printf "%-30s: %-8s : %s\\n" "AMRFINDER" "FAILED" "${amr_file} does not exist"  >> "${sample_name}.synopsis"
  status="FAILED"
fi

if [[ "${internal_phoenix}" == "true" ]]; then
  if [[ "${run_type}" == "all" ]]; then
    # check SRST2 output
    if [[ -s "${srst2_file}" ]]; then
      srst2_DB=$(echo "${srst2_file}" | rev | cut -d'_' -f4,5 | rev)
      info_NARC_List=$(head -n 1 "${srst2_file}")
      IFS='	' read -r -a NARC_array <<< "${info_NARC_List}"
      NARC_Num="${#NARC_array[@]}"
      NARC_Num=$(( NARC_Num - 1 ))
      #echo "${info_ResGANNCBI_List} - ${ResGANNCBI_Num}"
      if [[ "${NARC_Num}" -eq 0 ]]; then
        printf "%-30s: %-8s : %s\\n" "SRST2" "ALERT" "Completed, but NO KNOWN AMR genes present from ${srst2_DB}"  >> "${sample_name}.synopsis"
        if [[ "${status}" == "SUCCESS" ]]; then
          status="ALERT"
        fi
      else
        printf "%-30s: %-8s : %s\\n" "SRST2" "SUCCESS" "${NARC_Num} gene(s) found from ${srst2_DB}"  >> "${sample_name}.synopsis"
      fi
    else
      printf "%-30s: %-8s : %s\\n" "SRST2" "FAILED" "${sample_name}__fullgenes__*.txt file does not exist"  >> "${sample_name}.synopsis"
    fi
  else
    printf "%-30s: %-8s : %s\\n" "SRST2" "NA" "Assembly only isolate"  >> "${sample_name}.synopsis"
  fi
fi

# check Replicons
Rep_genes="NA"
if [[ -s "${gamma_replicon}" ]]; then
  Rep_genes=0
  repDB=$(echo ${gamma_replicon} | rev | cut -d'_' -f2,3 | rev)
  while read line_in; do
    if [[ "${line_in}" = "Gene	Contig	Start	Stop	Match_Type"* ]]; then
      :
    else
      Rep_genes=$(( Rep_genes + 1 ))
    fi
  done < "${gamma_replicon}"
  if [[ ${Rep_genes} -eq 0 ]]; then
    printf "%-30s: %-8s : %s\\n" "PLASMID_REPLICONS" "SUCCESS" "No replicons were found from ${repDB}"  >> "${sample_name}.synopsis"
  elif [[ ${Rep_genes} -ge 1 ]]; then
    printf "%-30s: %-8s : %s\\n" "PLASMID_REPLICONS" "SUCCESS" "${Rep_genes} replicon(s) found from ${repDB}"  >> "${sample_name}.synopsis"
  else
    echo "Should never get here (Rep_genes less than 0)"
  fi
else
  printf "%-30s: %-8s : %s\\n" "PLASMID_REPLICONS" "FAILED" "${sample_name}_PF-*.gamma does not exist"  >> "${sample_name}.synopsis"
  status="FAILED"
fi

# check HV
HV_genes="NA"
if [[ -s "${gamma_HV}" ]]; then
  HV_genes=0
  HVDB=$(echo ${gamma_HV} | rev | cut -d'_' -f2,3 | rev)
  while read line_in; do
    if [[ "${line_in}" = "Gene	Contig	Start	Stop	Match_Type"* ]]; then
      :
    else
      HV_genes=$(( HV_genes + 1 ))
    fi
  done < "${gamma_HV}"
  if [[ ${HV_genes} -eq 0 ]]; then
    printf "%-30s: %-8s : %s\\n" "HYPERVIRULENCE" "SUCCESS" "No hypervirulence genes were found from ${HVDB}"  >> "${sample_name}.synopsis"
  elif [[ ${HV_genes} -ge 1 ]]; then
    printf "%-30s: %-8s : %s\\n" "HYPERVIRULENCE" "SUCCESS" "${HV_genes} gene(s) found from ${HVDB}"  >> "${sample_name}.synopsis"
  else
    echo "Should never get here (HV_genes less than 0)"  >> "${sample_name}.synopsis"  >> "${sample_name}.synopsis"
  fi
else
  printf "%-30s: %-8s : %s\\n" "HYPERVIRULENCE" "FAILED" "${sample_name}_Hyper_Virulence_*.gamma does not exist"  >> "${sample_name}.synopsis"
  status="FAILED"
fi

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

if [[ "${internal_phoenix}" == "true" ]]; then
  printf "\n*BUSCO defines core genes as single-copy orthologs that should be highly conserved among the closely related species."  >> "${sample_name}.synopsis"
fi



#Script exited gracefully (unless something else inside failed)
exit 0
