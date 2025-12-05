process SCAFFOLDS_SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    stageInMode 'copy'
    // base_v2.3.0 - MUST manually change below (line 21)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    path samplesheet

    output:
    path '*.valid.csv' , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def container_version = "base_v2.3.0"
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
