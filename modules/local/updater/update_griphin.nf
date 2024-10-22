process UPDATE_GRIPHIN {
    tag "${project_id}"
    label 'process_low'
    stageInMode 'copy' // you need this or openpyxl complains that excel files aren't excel files. 
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    path(griphins_excel)
    //path(griphins_tsv)
    path(outdir) // output directory used as prefix for the summary file
    val(project_id)
    val(coverage)

    output:
    path("${project_id}_GRiPHin_Summary.xlsx"),    emit: griphin_report
    path("${project_id}_GRiPHin_Summary.tsv"),     emit: griphin_tsv_report
    path("versions.yml"),              emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    if (griphins_excel.size() == 2) {
        // Case where only two files are passed
        griphin_input = "-g1 ${griphins_excel[0]} -g2 ${griphins_excel[1]}"
    } else {
        // Case where griphins_excel contains many
        griphin_input = "--griphin_list"
    }
    """
    ${ica}combine_GRiPHins.py ${griphin_input} --output ${project_id}_GRiPHin_Summary --coverage ${coverage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
       combine_GRiPHins.py: \$(${ica}combine_GRiPHins.py --version)
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}