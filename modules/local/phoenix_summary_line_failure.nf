process CREATE_SUMMARY_LINE_FAILURE {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(synopsis), \
    path(fastp_total_qc), \
    path(trimd_ksummary), \
    path(taxonomy_file), \
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
        -x $taxonomy_file \\
        -o ${prefix}_summaryline.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
