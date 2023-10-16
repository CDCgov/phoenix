process GATHER_SUMMARY_LINES {
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path(summary_line_files)
    path(outdir_path)
    val(busco_val)

    output:
    path('Phoenix_Summary.tsv'), emit: summary_report
    path("versions.yml")             , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def busco_parameter = busco_val ? "--busco" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}Create_phoenix_summary_tsv.py \\
        --out Phoenix_Summary.tsv \\
        $busco_parameter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
