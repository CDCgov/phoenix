process SRATOOLS_PREFETCH {
    tag "$id"
    label 'process_low'

    conda (params.enable_conda ? 'bioconda::sra-tools=2.11.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5321ha49a11a_3' :
        'quay.io/biocontainers/sra-tools:2.11.0--pl5321ha49a11a_3' }"

    input:
    path(id)
    

    output:
    path(id)                 , emit: sra
    path 'versions.yml'      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    prefetch --option-file $id --type fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}
