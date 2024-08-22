process UPDATE_GRIPHIN {
    tag "${meta.project_id}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.0.2'

    input:
    tuple val(meta), path(old_griphin), path(new_griphin), path(outdir) // output directory used as prefix for the summary file
    val(coverage)

    output:
    path("*_GRiPHin_Summary.xlsx"),    emit: griphin_report
    path("*_GRiPHin_Summary.tsv"),     emit: griphin_tsv_report
    path("versions.yml"),              emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}combine_GRiPHins.py -g1 ${old_griphin} -g2 ${new_griphin} --output ${outdir}_GRiPHin_Summary --coverage ${coverage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
       combine_GRiPHins.py: \$(${ica}combine_GRiPHins.py --version)
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}