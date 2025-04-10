process CIRCLATOR {
    tag "${meta.id}"
    label 'process_medium'
    container 'staphb/circlator:1.5.6'
    //sha256:04576d9de1dae0244e96e69f4bc8f1f943849e6967733ce281cbcdaa5706f16f
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_unpolished.fasta"), emit: fasta
    //tuple val(meta), path("*_unpolished.log"),   emit: circlator_log
    path ("versions.yml"),                       emit: versions

    script:
    """
    # is unpolished the right name for this?
    circlator fixstart $fasta ${meta.id}_unpolished

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circlator: \$( circlator --version | sed -e "s/CIRCLATOR v//g" )
    END_VERSIONS
    """
}