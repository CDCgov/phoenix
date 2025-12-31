process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    stageInMode 'copy'
    // base_v2.2.0 - MUST manually change below (line 24)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    path(samplesheet)
    val(reads_entry)
    val(scaffolds_entry)
    val(directory_entry)
    val(meta) // used for --mode update_phoenix to get meta.full_project_id - to make sure things are published to the right dir in --input

    output:
    path('samplesheet.valid.csv'),                   emit: csv
    path('samplesheet.valid_*.csv'),  optional:true, emit: csv_by_dir // only need if multiple dirs in --input for --pipeline update_phoenix
    path("versions.yml"),                            emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def container_version = "base_v2.2.0"
    def reads_check = reads_entry ? "true" : "false"
    def scaffolds_check = scaffolds_entry ? "true" : "false"
    def directory_check = directory_entry ? "true" : "false"
    def sheet_by_dir = (params.mode_upper == "UPDATE_PHOENIX" || params.mode_upper == "CENTAR") ? "--sheet_by_dir" : ""
    def updater = (params.mode_upper == "UPDATE_PHOENIX" || params.mode_upper == "CENTAR") ? "--updater" : ""
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
        ${ica}check_directory_samplesheet.py ${samplesheet} samplesheet.valid.csv ${updater} ${sheet_by_dir}
        script_version=\$(echo check_directory_samplesheet.py: \$(${ica}check_directory_samplesheet.py --version ))
    else
        echo "No valid check type provided, exiting."
        exit 1
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