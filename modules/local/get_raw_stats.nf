process GET_RAW_STATS {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(reads)


    output:
    tuple val(meta), path('*_stats.txt'),optional:true,                     emit: raw_stats
    tuple val(meta), path('*_raw_read_counts.txt'),optional:true,           emit: combined_raw_stats
    path("versions.yml"),optional:true,                                     emit: versions
    tuple val(meta), path('*_result.txt'),optional:true,                    emit: outcome
    tuple val(meta), path('*.fastq.gz'),optional:true,                      emit: reads
    tuple val(meta),path('*_prdresult.txt'),optional:true,                  emit: compr_rds
    tuple val(meta),path('*_summaryline_failure.tsv'),optional:true,        emit: summary_line
    
    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num1 = "${reads[0]}".minus(".fastq.gz")
    def num2 = "${reads[1]}".minus(".fastq.gz")

    """
    set +e

    # check for file corruption

    for isolate in ${reads}
        do
            gzip -t ${reads[0]} 2>> ${num1}.txt
            gzip -t ${reads[1]} 2>> ${num2}.txt
        done

    if grep "error" ${num1}.txt || grep "error" ${num2}.txt || grep "unexpected" ${num1}.txt || grep "unexpected" ${num2}.txt; then
        echo "FAILED CORRUPTION CHECK! CANNOT UNZIP FASTQ FILE. CHECK FASTQ FILE(S) FOR CORRUPTION!" > ${prefix}_results.txt
        # error warning for line_summary channel
        echo "ID	Auto_QC_Outcome	Warning_Count	Estimated_Coverage	Genome_Length	Assembly_Ratio_(STDev)	#_of_Scaffolds_>500bp	GC_%	Species	Taxa_Confidence	Taxa_Coverage	Taxa_Source	Kraken2_Trimd	Kraken2_Weighted	MLST_Scheme_1	MLST_1	MLST_Scheme_2	MLST_2	GAMMA_Beta_Lactam_Resistance_Genes	GAMMA_Other_AR_Genes	AMRFinder_Point_Mutations	Hypervirulence_Genes	Plasmid_Incompatibility_Replicons	Auto_QC_Failure_Reason" > ${prefix}_summaryline_failure.tsv
        #file contents
        echo "${prefix}	FAIL	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	kraken2_trimmed	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown		CANNOT UNZIP FASTQ FILE. CHECK FASTQ FILE(S) FOR CORRUPTION!" | tr -d '\n' >> ${prefix}_summaryline_failure.tsv
    else
        echo "PASS" > ${prefix}_results.txt
    fi
    
    # proceed to cumulative read counts if files aren't corrupt
    if grep -Fx "PASS" ${prefix}_results.txt
        then
        q30.py ${reads[0]} > ${prefix}_R1_stats.txt
        q30.py ${reads[1]} > ${prefix}_R2_stats.txt
    fi

    if [ -f ${prefix}_R2_stats.txt -a -f ${prefix}_R1_stats.txt ] 
        then
        create_raw_stats_output.py -n ${prefix} -r1 ${prefix}_R1_stats.txt -r2 ${prefix}_R2_stats.txt
        comb_stats_chk.py -r ${prefix}_raw_read_counts.txt
    fi

    if grep "PASS" ${prefix}_result.txt
    then
        mv ${reads[0]} ${num1}_C.fastq.gz
        mv ${reads[1]} ${num2}_C.fastq.gz
    else
        if [ ! -f ${prefix}_summaryline_failure.tsv ]
        then
            # error warning for line_summary channel
            echo "ID	Auto_QC_Outcome	Warning_Count	Estimated_Coverage	Genome_Length	Assembly_Ratio_(STDev)	#_of_Scaffolds_>500bp	GC_%	Species	Taxa_Confidence	Taxa_Coverage	Taxa_Source	Kraken2_Trimd	Kraken2_Weighted	MLST_Scheme_1	MLST_1	MLST_Scheme_2	MLST_2	GAMMA_Beta_Lactam_Resistance_Genes	GAMMA_Other_AR_Genes	AMRFinder_Point_Mutations	Hypervirulence_Genes	Plasmid_Incompatibility_Replicons	Auto_QC_Failure_Reason" > ${prefix}_summaryline_failure.tsv
            #file contents
            echo "${prefix}	FAIL	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	kraken2_trimmed	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown		Read pairs are NOT the same!" | tr -d '\n' >> ${prefix}_summaryline_failure.tsv
        fi
        # delete files that would be enter the channel
        rm ${reads[0]}
        rm ${reads[1]}
        mv ${prefix}_result.txt ${prefix}_prdresult.txt
        rm ${prefix}_prdresult.txt

    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}