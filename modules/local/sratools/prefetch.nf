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
    path 'tmp_inputs/*'             , emit: sra
    path 'sra_samples.csv'          , emit: samplesheet
    path 'versions.yml'      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    TMP_INPUTS=tmp_inputs
    mkdir "\$TMP_INPUTS"

    # fetch sras
    prefetch --option-file $id --output-directory "\$TMP_INPUTS"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}