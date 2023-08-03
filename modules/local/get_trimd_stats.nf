process GET_TRIMD_STATS {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.0.2'

    input:
    tuple val(meta), path(fastp_trimd_json), path(fastp_singles_json)

    output:
    tuple val(meta), path('*_trimmed_read_counts.txt'), emit: fastp_total_qc
    path("versions.yml"),                               emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    FastP_QC.py \\
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
