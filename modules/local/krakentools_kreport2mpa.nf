def VERSION = '1.2' // Version information not provided by tool on CLI

process KRAKENTOOLS_KREPORT2MPA {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(kraken_report)

    output:
    tuple val(meta), path('*.mpa'), emit: mpa
    path("versions.yml")          , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}kreport2mpa.py \\
        --report-file ${kraken_report} \\
        --output ${prefix}.mpa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools: $VERSION
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}