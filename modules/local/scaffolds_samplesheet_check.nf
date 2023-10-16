process SCAFFOLDS_SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path samplesheet

    output:
    path '*.valid.csv' , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}check_assembly_samplesheet.py \\
    $samplesheet \\
    samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container} 
    END_VERSIONS
    """
}
