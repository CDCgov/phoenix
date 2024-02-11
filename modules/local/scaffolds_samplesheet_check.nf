process SCAFFOLDS_SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.0.2'

    input:
    path samplesheet

    output:
    path '*.valid.csv' , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    check_assembly_samplesheet.py \\
    $samplesheet \\
    samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container} 
    END_VERSIONS
    """
}
