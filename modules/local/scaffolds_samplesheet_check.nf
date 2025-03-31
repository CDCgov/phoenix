process SCAFFOLDS_SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    stageInMode 'copy'
    // base_v2.2.0 - MUST manually change below (line 21)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    path samplesheet

    output:
    path '*.valid.csv' , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}check_assembly_samplesheet.py \\
    $samplesheet \\
    samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        check_assembly_samplesheet.py: \$(${ica}check_assembly_samplesheet.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container} 
    END_VERSIONS
    """
}
