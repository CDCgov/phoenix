process READ_COUNT_CHECK {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(raw_read_counts), path(fairy_corrupt_outcome)
    val(busco_val)

    output:
    tuple val(meta), path('*_results.txt'),                emit: outcome 
    path('*_summaryline_failure.tsv'),      optional:true, emit: summary_line
    tuple val(meta), path('*_results_old_2.txt'),          emit: outcome_to_edit
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
    def busco_parameter = busco_val ? "--busco" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    # Output check for messages indicating read pairs that do not match
    ${ica}fairy.py -r ${raw_read_counts} -f ${fairy_corrupt_outcome} ${busco_parameter}

    #making a copy of the results file to pass to BBMAP_REFORMAT to handle file names being the same
    cp ${prefix}_results.txt ${prefix}_results_old_2.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}