process FILE_RENAME {
    label 'process_single'
    container parmas.phoenix_base_container

    input:
    path(griphins)

    output:
    path('*_GRiPHin.*'),  emit: renamed_griphins
    path("versions.yml"), emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def container_version = params.phoenix_container_version
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}file_rename.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
        \$(${ica}file_rename.py --version)
    END_VERSIONS
    """
}