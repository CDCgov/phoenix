process GATHER_SUMMARY_LINES {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    path(summary_line_files)

    output:
    path 'Phoenix_Output_Report.tsv'  , emit: summary_report
    path "versions.yml"                , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """
    Create_phoenix_summary_tsv.py \\
        --out Phoenix_Output_Report.tsv \\
        $summary_line_files

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
