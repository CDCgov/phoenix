process CREATE_SUMMARY_LINE_FAILURE {
    tag "${meta.id}"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(synopsis), \
    path(fastp_total_qc), \
    path(trimd_ksummary), \
    val(spades_outcome)
    
    output:
    path('*_summaryline.tsv') , emit: line_summary
    path "versions.yml"       , emit: versions

    when:
    "${spades_outcome[0]}" == "run_failure" || "${spades_outcome[1]}" == "no_scaffolds" || "${spades_outcome[2]}" == "no_contigs"

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    Phoenix_summary_line.py \\
        -n ${prefix} \\
        -k $trimd_ksummary \\
        -t $fastp_total_qc \\
        -s $synopsis \\
        -o ${prefix}_summaryline.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
