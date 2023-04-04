def VERSION = '1.2' // Version information not provided by tool on CLI

process KRAKEN2_KRONA {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

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