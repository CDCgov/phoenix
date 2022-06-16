process CREATE_SUMMARY_LINE {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(trimmed_qc_data_file)
    tuple val(meta), path(mlst_file)
    tuple val(meta), path(hypervirulence_gamma_file)
    tuple val(meta), path(ar_gamma_file)
    tuple val(meta), path(quast_report)
    tuple val(meta), path(ratio_file)

    output:
    path '*_summaryline.tsv'           , emit: line_summary
    path "versions.yml"                , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Phoenix_summary_line.py \\
        -q $quast_report \\
        -t $trimmed_qc_data_file \\
        -a $ar_gamma_file \\
        -v $hypervirulence_gamma_file \\
        -r $ratio_file \\
        -m $mlst_file \\
        -n ${prefix} \\
        -o ${prefix}_summaryline.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
