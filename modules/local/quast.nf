process QUAST {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::quast=5.0.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2' :
        'quay.io/biocontainers/quast:5.0.2--py37pl526hb5aa323_2' }"

    input:
    tuple val(meta), path(consensus)

    output:
    tuple val(meta), path("quast")        , emit: results
    tuple val(meta), path('*.tsv')        , emit: report_tsv
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    quast.py \\
        --output-dir quast \\
        --threads $task.cpus \\
        $args \\
        $consensus

    mv quast/report.tsv ./${prefix}_report.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
