process GATHERING_READ_QC_STATS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/python:3.8.3' :
      'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(fastp_trimd_json), path(fastp_singles_json)

    output:
    tuple val(meta), path('*_raw_read_counts.txt')     , emit: fastp_raw_qc
    tuple val(meta), path('*_trimmed_read_counts.txt') , emit: fastp_total_qc

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    FastP_QC.py \\
      --trimmed_json ${fastp_trimd_json} \\
      --single_json ${fastp_singles_json} \\
      --name ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
