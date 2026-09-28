process GATHER_SUMMARY_LINES {
    label 'process_single'
    container params.phoenix_base_container

    input:
    val(meta) // need for meta.full_project_id in -profile update_phoenix /species specific pipeliens for publishing
    path(summary_line_files)
    path(full_project_id) // for --mode update_phoenix and species specific pipelines for publishing location
    val(busco_val)
    val(pipeline_info) //needed for --pipeline update_phoenix

    output:
    path('*hoenix_Summary.tsv'), emit: summary_report
    path("versions.yml")       , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def busco_parameter = busco_val ? "--busco" : ""
    def software_versions = pipeline_info ? "--software_versions ${pipeline_info}" : ""
    def container_version = params.phoenix_container_version
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def updater = ((params.mode_upper == "UPDATE_PHOENIX" && params.outdir == "${launchDir}/phx_output") || (params.mode_upper == "CENTAR" && params.outdir == "${launchDir}/phx_output")) ? "--all_samples" : ""
    def output = "Phoenix_Summary.tsv" 
    """
    ${ica}Create_phoenix_summary_tsv.py --out ${output} ${updater} ${software_versions} $busco_parameter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        \$(${ica}Create_phoenix_summary_tsv.py --version)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}