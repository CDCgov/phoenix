process GATHER_SUMMARY_LINES {
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 22)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    val(meta) // need for meta.full_project_id in -profile update_phoenix /species specific pipeliens for publishing
    path(summary_line_files)
    path(outdir_path)
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
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def updater = ((params.mode_upper == "UPDATE_PHOENIX" && params.outdir == "${launchDir}/phx_output") || (params.mode_upper == "CENTAR" && params.outdir == "${launchDir}/phx_output")) ? "--all_samples" : ""
    def output = "Phoenix_Summary.tsv" 
    """
    ${ica}Create_phoenix_summary_tsv.py --out ${output} ${updater} ${software_versions} $busco_parameter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        Create_phoenix_summary_tsv.py: \$(${ica}Create_phoenix_summary_tsv.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
