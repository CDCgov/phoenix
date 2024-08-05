process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    stageInMode 'copy'
    // base_v2.1.0 - MUST manually change below (line 20)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    path(samplesheet)
    val(reads_entry)
    val(scaffolds_entry)
    val(directory_entry)

    output:
    path('samplesheet.valid.csv'),  emit: csv
    path("versions.yml"),           emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container_version = "base_v2.1.0"
    def reads_check = reads_entry ? "true" : "false"
    def scaffolds_check = scaffolds_entry ? "true" : "false"
    def directory_check = directory_entry ? "true" : "false"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    if [ ${reads_check} = "true" ]; then
        echo "Running check of reads samplesheet"
        ${ica}check_samplesheet.py ${samplesheet} samplesheet.valid.csv
        script_version=\$(echo check_samplesheet.py: \$(${ica}check_samplesheet.py --version ))
    elif [ ${scaffolds_check} = "true" ]; then
        echo "Running check of assembly samplesheet"
        ${ica}check_assembly_samplesheet.py ${samplesheet} samplesheet.valid.csv
        script_version=\$(echo check_assembly_samplesheet.py: \$(${ica}check_assembly_samplesheet.py --version ))
    elif [ ${directory_check} = "true" ]; then
        echo "Running check of directory samplesheet"
        ${ica}check_directory_samplesheet_v2.py ${samplesheet} samplesheet.valid.csv
        script_version=\$(echo check_directory_samplesheet_v2.py: \$(${ica}check_directory_samplesheet.py --version ))
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        \${script_version}
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container} 
    END_VERSIONS
    """
}