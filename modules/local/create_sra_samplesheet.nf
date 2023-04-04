process CREATE_SRA_SAMPLESHEET {
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(samplesheet) // 
    path(metadata_csv)
    path(directory)

    output:
    path('samplesheet.csv'), emit: csv
    path("versions.yml"),    emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """
    sra_samplesheet.py -d ${directory}/raw_fastqs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}