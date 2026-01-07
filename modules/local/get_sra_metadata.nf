process ENTREZDIRECT_ESEARCH {
    tag "${meta.id}"
    label 'process_single'
    maxForks 3
    // v24.0--he881be0_0
    container 'quay.io/biocontainers/entrez-direct@sha256:ccafbda537a8ab77206758c91a383defe0ea5007365b526aa89db3cfbf451d51'

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