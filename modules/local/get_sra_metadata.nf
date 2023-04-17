process ENTREZDIRECT_ESEARCH {
    tag "${sra_folder}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/entrez-direct:16.2--he881be0_1':
        'quay.io/biocontainers/entrez-direct:16.2--he881be0_1' }"

    input:
    path(sra_folder) // using the folder name to get the srr number
    val(database)

    output:
    path("*_sra_metadata.csv"), emit: metadata_csv
    path("versions.yml"),       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def srr_number = sra_folder.toString()
    """
    esearch \\
        -db $database \\
        -query $srr_number | \\
        efetch -format runinfo > ${srr_number}_sra_metadata.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        esearch: \$(esearch -version)
    END_VERSIONS
    """
}