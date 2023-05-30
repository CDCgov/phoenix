process GET_RAW_STATS {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(reads)


    output:
    tuple val(meta), path('*_stats.txt'),                          emit: raw_stats
    tuple val(meta), path('*_raw_read_counts.txt'),                emit: combined_raw_stats
    path("versions.yml"),                                          emit: versions
    tuple val(meta), path('*_result.txt'),                         emit: outcome
    tuple val(meta), path('*.fastq.gz'),optional:true,             emit: reads

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num1 = "${reads[0]}".minus(".fastq.gz")
    def num2 = "${reads[1]}".minus(".fastq.gz")

    """
    # check for file corruption
    do
        gzip -t ${reads[0]} 2>> ${num1}.txt
        gzip -t ${reads[1]} 2>> ${num2}.txt
    done

    # may be able to bypass
    if grep "error" ${num1}.txt || grep "error" ${num2}.txt || grep "unexpected" ${num1}.txt || grep "unexpected" ${num2}.txt; then
        echo "FAILED CORRUPTION CHECK! CANNOT UNZIP FASTQ FILE. CHECK FASTQ FILE(S) FOR CORRUPTION!" > ${prefix}_results.txt
    else
        echo "PASS" > ${prefix}_results.txt
    fi
    
    # proceed to cumulative read counts if files aren't corrupt
    if grep -Fx "PASS" ${prefix}_results.txt
        then
        q30.py ${reads[0]} > ${prefix}_R1_stats.txt
        q30.py ${reads[1]} > ${prefix}_R2_stats.txt
    else 
        mv ${prefix}_results.txt ${prefix}_raw_read_counts.txt
    fi

    if [ -f ${prefix}_R2_stats.txt -a -f ${prefix}_R1_stats.txt ] 
        then
        create_raw_stats_output.py -n ${prefix} -r1 ${prefix}_R1_stats.txt -r2 ${prefix}_R2_stats.txt
        comb_stats_chk.py -r ${prefix}_raw_read_counts.txt
    else echo "YOUR READ PAIRS ARE NOT THE SAME! THESE SAMPLES HAVE BEEN SKIPPED. PHOENIX ONLY ANALYZES ISOLATES WITH THE SAME NUMBER OF READS!" > ${prefix}_raw_read_counts.txt
    fi

    # only send the reads that pass all QC checks
    if grep -Fx "PASS" ${prefix}_result.txt
        then
        mv ${reads[0]} ${num1}_C.fastq.gz
        mv ${reads[1]} ${num2}_C.fastq.gz
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}