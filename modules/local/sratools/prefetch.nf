process SRATOOLS_PREFETCH {
    tag "$sra_samples"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5321ha49a11a_3' :
        'quay.io/biocontainers/sra-tools:2.11.0--pl5321ha49a11a_3' }"

    input:
    path(sra_samples)

    output:
    path('SRR*')             , emit: sra_folder
    path('sra_samples.csv')  , emit: samplesheet
    path('versions.yml')     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # fetch sras
    prefetch --option-file $sra_samples

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}