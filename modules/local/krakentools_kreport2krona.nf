def VERSION = 'https://github.com/jenniferlu717/KrakenTools/commit/1271ae2ee2289148f9d4bae4a59323d7a8ea288a' // Version information not provided by tool on CLI

process KRAKEN2_KRONA {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.0.0'

    /*container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"*/

    input:
    tuple val(meta), path(kraken_report)
    val(type) //weighted, trimmmed or assembled

    output:
    tuple val(meta), path('*.krona'), emit: krona
    path "versions.yml"             , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/ orginally from https://github.com/jenniferlu717/KrakenTools on 6/15/2022
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kreport2krona.py \\
        --report ${kraken_report} \\
        --output ${prefix}_${type}.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools: $VERSION
    END_VERSIONS
    """
}