def VERSION = '1.2' // Version information not provided by tool on CLI

process KRAKENTOOLS_KREPORT2MPA {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(kraken_report)

    output:
    tuple val(meta), path('*.mpa'), emit: mpa
    path "versions.yml"           , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kreport2mpa.py \\
        --report-file ${kraken_report} \\
        --output ${prefix}.mpa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools: $VERSION
    END_VERSIONS
    """
}