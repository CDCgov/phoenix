process CREATE_SAMPLESHEET {
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.0.2'

    input:
    path(directory)

    output:
    path("GRiPHin_samplesheet_created.csv"), emit: samplesheet
    path("versions.yml"),                    emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/griphin/bin/
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    create_samplesheet.py --directory ${directory}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}