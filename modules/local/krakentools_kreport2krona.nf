def VERSION = '1.2' // Version information not provided by tool on CLI

process KRAKEN2_KRONA {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(kraken_report)
    val(type) //weighted, trimmmed or assembled

    output:
    tuple val(meta), path('*.krona'), emit: krona
    path("versions.yml")            , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/ orginally from https://github.com/jenniferlu717/KrakenTools on 6/15/2022
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) {
        ica = ""
    } else if (params.ica==true) {
        ica = "python ${workflow.launchDir}/bin/"
    } else {
        error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}kreport2krona.py \\
        --report ${kraken_report} \\
        --output ${prefix}_${type}.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools: $VERSION
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}