process NANOQ {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/nanoq:0.10.0--h031d066_2'
    //sha256:e3f7fc6e04ed0b2ae8753264c9898d981f798ada6a41689bf788e40824816ae4

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_trim.fastq.gz"), emit: fastq
    tuple val(meta), path("*_stats_mqc.csv"), emit: nano_stats
    path ("versions.yml"),                    emit: versions

    script:
    """
    nanoq -i $reads --min-len 1000 --min-qual 15 -r ${meta.id}_stats_mqc.txt --stats --header -o ${meta.id}_trim.fastq.gz >> ${meta.id}_nanoq.log

    # convert to stats for reporting in griphin later to csv
    sed 's/ /,/g' ${meta.id}_stats_mqc.txt > ${meta.id}_stats_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e 's/nanoq //g')
    END_VERSIONS
    """

}
