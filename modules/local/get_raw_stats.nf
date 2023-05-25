process GET_RAW_STATS {
    tag "$meta.id"
    label 'process_single'
    label 'error_ignore'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(reads)


    output:
    tuple val(meta), path('*_stats.txt'),           emit: raw_stats
    tuple val(meta), path('*_raw_read_counts.txt'), emit: combined_raw_stats
    path("versions.yml"),                           emit: versions
    tuple val(meta), path('*_result.txt'),          emit: outcome
    tuple val(meta), path('*.fastq.gz'),             emit: reads

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num1 = "${reads[0]}".minus(".fastq.gz")
    def num2 = "${reads[1]}".minus(".fastq.gz")

    """
    q30.py ${reads[0]} > ${prefix}_R1_stats.txt
    q30.py ${reads[1]} > ${prefix}_R2_stats.txt

    
    if [ -f ${prefix}_R2_stats.txt -a -f ${prefix}_R1_stats.txt ]; then
        create_raw_stats_output.py -n ${prefix} -r1 ${prefix}_R1_stats.txt -r2 ${prefix}_R2_stats.txt
        comb_stats_chk.py -r ${prefix}_raw_read_counts.txt

        mv ${reads[0]} ${num1}C.fastq.gz
        mv ${reads[1]} ${num2}C.fastq.gz
    
        mv ${num1}C.fastq.gz ${num1}.fastq.gz
        mv ${num2}C.fastq.gz ${num2}.fastq.gz
        else echo "FAIL" > ${prefix}_raw_read_counts.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}