process NANOQ {
    tag "$meta"
    label 'process_medium'
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/nanoq:0.10.0--h031d066_2'
    //sha256:e3f7fc6e04ed0b2ae8753264c9898d981f798ada6a41689bf788e40824816ae4
    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz") , emit: fastq
    tuple val(meta), path("*_stats_mqc.csv") , emit: nano_summary
    path "versions.yml"                       , emit: versions

    script:
    """
    nanoq -i $reads -l 1000 -q 15 -r ${meta}_stats_mqc.txt -s -H -o ${meta}_trim.fastq.gz >> nanoq.log 
    sed 's/ /,/g' ${meta}_stats_mqc.txt > ${meta}_stats_mqc.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$( nanoq --version | sed -e "s/NANOQ v//g" )
    END_VERSIONS

    """

}
