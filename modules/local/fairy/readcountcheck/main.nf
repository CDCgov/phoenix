process READ_COUNT_CHECK {
    tag "${meta.id}"
    label 'process_medium'
    // base_v2.1.0 - MUST manually change below (line 26)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(raw_read_counts), path(fairy_corrupt_outcome)
    val(busco_val)

    output:
    tuple val(meta), path('*_summary.txt'),                emit: outcome
    path('*_summaryline.tsv'),      optional:true, emit: summary_line
    tuple val(meta), path('*_summary_old_2.txt'),          emit: outcome_to_edit
    tuple val(meta), path('*.synopsis'),    optional:true, emit: synopsis
    path("versions.yml"),                                  emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" }
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_parameter = busco_val ? "--busco" : ""
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    # Output check for messages indicating read pairs that do not match
    ${ica}fairy.py -r ${raw_read_counts} -f ${fairy_corrupt_outcome} ${busco_parameter}

    # adding an extra line to allow KRAKEN2_TRIMD work with the when statement. Will overwrite in SCAFFOLDS_CHECK_COUNT module
    echo -e "\\nEnd_of_File" >> ${prefix}_summary.txt

    #making a copy of the summary file to pass to BBMAP_REFORMAT to handle file names being the same
    cp ${prefix}_summary.txt ${prefix}_summary_old_2.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        fairy.py: \$(${ica}fairy.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
