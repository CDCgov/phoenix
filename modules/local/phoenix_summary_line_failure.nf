process CREATE_SUMMARY_LINE_FAILURE {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 30)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

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
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extended_qc_arg = extended_qc ? "--extended_qc" : ""
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script = params.ica ? "python ${params.ica_path}/Phoenix_summary_line.py" : "Phoenix_summary_line.py"
    """
    ${script} \\
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
        Phoenix_summary_line.py: \$(${script} --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
