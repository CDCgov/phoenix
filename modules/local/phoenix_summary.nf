process GATHER_SUMMARY_LINES {
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 22)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    path(summary_line_files)
    path(outdir_path)
    val(busco_val)

    output:
    path('Phoenix_Summary.tsv'), emit: summary_report
    path("versions.yml")       , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def busco_parameter = busco_val ? "--busco" : ""
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Create_phoenix_summary_tsv.py \\
        --out Phoenix_Summary.tsv \\
        $busco_parameter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        Create_phoenix_summary_tsv.py: \$(${ica}Create_phoenix_summary_tsv.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
