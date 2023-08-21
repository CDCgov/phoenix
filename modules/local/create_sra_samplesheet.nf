process CREATE_SRA_SAMPLESHEET {
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path(renamed_reads)
    path(metadata_csvs)
    path(directory)
    val(srr_param)

    output:
    path('sra_samplesheet.csv'), emit: csv
    path("versions.yml"),        emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) {
        ica = ""
    } else if (params.ica==true) {
        ica = "python ${workflow.launchDir}/bin/"
    } else {
        error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods."
    }
    def use_srr = srr_param ? "--use_srr" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    full_path=\$(readlink -f ${directory})

    ${ica}sra_samplesheet.py -d \$full_path $use_srr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}