process ENTREZDIRECT_ESEARCH {
    tag "${meta.id}"
    label 'process_single'
    maxForks 3
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    tuple val(meta), path(sra_folder)

    output:
    tuple val(meta), path("*_sra_metadata.csv"), emit: metadata_csv
    path("versions.yml"),                        emit: versions

    script:
    """
    esearch \\
        -db sra \\
        -query ${meta.id} | \\
        efetch -format runinfo > ${meta.id}_sra_metadata.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    END_VERSIONS
    """
}