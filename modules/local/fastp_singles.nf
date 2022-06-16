process FASTP {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::fastp=0.23.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h79da9fb_0' :
        'quay.io/biocontainers/fastp:0.23.2--h79da9fb_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.singles.fastq.gz')  , emit: reads
    tuple val(meta), path('*.json')              , emit: json
    tuple val(meta), path('*.html')              , emit: html
    tuple val(meta), path('*.log')               , emit: log
    path "versions.yml"                          , emit: versions
    tuple val(meta), path('*.merged.fastq.gz')   , optional:true, emit: reads_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat ${reads[0]} ${reads[1]} > ${prefix}.cat_singles.fastq
    gzip ${prefix}.cat_singles.fastq
    fastp \\
        --in1 ${prefix}.cat_singles.fastq.gz \\
        --thread $task.cpus \\
        --json ${prefix}_singles.fastp.json \\
        --html ${prefix}_singles.fastp.html \\
        --out1 ${prefix}.singles.fastq.gz \\
        $args \\
        2> ${prefix}.fastp.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}