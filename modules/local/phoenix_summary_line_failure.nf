process CREATE_SUMMARY_LINE_FAILURE {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 30)!!!
    container 'quay.io/jvhagey/phoenix@sha256:caa2a5660c73d0376d7beb14069436a0e2403bda68904ff140cb789bf4f8753d'

    input:
    tuple val(meta), path(synopsis), \
    path(fastp_total_qc), \
    path(trimd_ksummary), \
    path(taxonomy_file), \
    val(spades_outcome)
    val(extended_qc)

    output:
    path('*_summaryline.tsv'), emit: line_summary
    path("versions.yml")     , emit: versions

    when:
    "${spades_outcome[0]}" == "run_failure" || "${spades_outcome[1]}" == "no_scaffolds" || "${spades_outcome[2]}" == "no_contigs"

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extended_qc_arg = extended_qc ? "--extended_qc" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Phoenix_summary_line.py \\
        -n ${prefix} \\
        -k $trimd_ksummary \\
        -t $fastp_total_qc \\
        -s $synopsis \\
        -x $taxonomy_file \\
        -o ${prefix}_summaryline.tsv \\
        $extended_qc_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        Phoenix_summary_line.py: \$(${ica}Phoenix_summary_line.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
