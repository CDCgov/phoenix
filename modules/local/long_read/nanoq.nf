process NANOQ {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/nanoq:0.10.0--h031d066_2'
    //sha256:e3f7fc6e04ed0b2ae8753264c9898d981f798ada6a41689bf788e40824816ae4

    input:
    tuple val(meta), path(rawstats)
    tuple val(meta), path(subfastq)

    output:
    tuple val(meta), path("*_trim.fastq.gz"),   emit: fastq
    tuple val(meta), path("*_nanoq_stats.csv"), emit: nano_stats
    tuple val(meta), path("*_nanoq.log"),       emit: nano_log
    path ("versions.yml"),                      emit: versions
    script:
    """
    nanoq -i $subfastq -l 2000 -q 15 -r trim.txt -s -H -o ${meta.id}_trim.fastq.gz >> ${meta.id}_nanoq.log 
    echo -e "raw_reads trim_reads bases n50 longest shortest mean_length median_length mean_quality median_quality\n\$(awk 'NR==2 {print \$4}' $rawstats) \$(awk 'NR==2 {print}' trim.txt)" > ${meta.id}_stats_mqc.txt
    # convert to stats for reporting in griphin later to csv
    sed 's/ /,/g' ${meta.id}_stats_mqc.txt > ${meta.id}_nanoq_stats.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e "s/nanoq //g")

    END_VERSIONS
    """

}
