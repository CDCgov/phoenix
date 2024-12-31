process ENTREZDIRECT_ESEARCH {
    tag "${meta.id}"
    label 'process_single'
    maxForks 3
    // v16.2--he881be0_1
    // container 'quay.io/biocontainers/entrez-direct@sha256:08a155e41ff29d396620a40906a16d3285fa21b508704ea161aeb3c2e071ef07'
    // v22.4--he881be0_0
    container 'quay.io/biocontainers/entrez-direct@sha256:07d6f38091457055aa9f93479fddb3a939eb703fc54445778521ae0b60966e6b'

    input:
    tuple val(meta), path(sra_folder)

    output:
    tuple val(meta), path("*_sra_metadata.csv"), emit: metadata_csv
    path("versions.yml"),                        emit: versions

    script:
    def container = task.container.toString() - "quay.io/biocontainers/entrez-direct@"
    """
    esearch \\
        -db sra \\
        -query ${meta.id} | \\
        efetch -format runinfo > ${meta.id}_sra_metadata.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
        esearch_container: ${container}
    END_VERSIONS
    """
}