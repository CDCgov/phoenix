process CIRCLATOR {
    tag "$meta"
    label 'process_medium'
    container 'staphb/circlator:1.5.6'
    //sha256:04576d9de1dae0244e96e69f4bc8f1f943849e6967733ce281cbcdaa5706f16f


    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fasta") , emit: fasta
    tuple val(meta), path("*mqc.tsv") , emit: circ_summary
    path "versions.yml"                       , emit: versions

    script:
    """
    circlator fixstart $fasta ${meta}_unpolished
    mv ${meta}_unpolished.log ${meta}_summary_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circlator: \$( circlator --version | sed -e "s/CIRCLATOR v//g" )
    END_VERSIONS
    """
}