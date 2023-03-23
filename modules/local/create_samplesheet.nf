process CREATE_SAMPLESHEET {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(directory)

    output:
    path("GRiPHin_samplesheet_created.csv"), emit: samplesheet
    path("versions.yml"),                    emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/griphin/bin/
    """
    create_samplesheet.py --directory ${directory}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}