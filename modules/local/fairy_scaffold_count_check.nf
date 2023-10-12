process SCAFFOLD_COUNT_CHECK {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*_results.txt'),                emit: outcome
    path('*_summaryline_failure.tsv'),      optional:true, emit: summary_line
    tuple val(meta), path('*.synopsis'),    optional:true, emit: synopsis
    path("versions.yml"),                                  emit: versions

    script:
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
    #checking that the output contains scaffolds still:
    if grep "Output:                 	0 reads (0.00%) 	0 bases (0.00%)" ${prefix}.bbmap_filtered.log; then
        echo "FAILED: No scaffolds left after filtering!" >> ${prefix}_result_Old.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
