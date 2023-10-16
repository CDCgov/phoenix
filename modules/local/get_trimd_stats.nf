process GET_TRIMD_STATS {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(fastp_trimd_json), path(fastp_singles_json)

    output:
    tuple val(meta), path('*_trimmed_read_counts.txt'), emit: fastp_total_qc
    path("versions.yml"),                               emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}FastP_QC.py \\
      --trimmed_json ${fastp_trimd_json} \\
      --single_json ${fastp_singles_json} \\
      --name ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
