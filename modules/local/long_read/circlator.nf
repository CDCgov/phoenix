process CIRCLATOR {
    tag "${meta.id}"
    label 'process_medium'
    container 'staphb/circlator@sha256:04576d9de1dae0244e96e69f4bc8f1f943849e6967733ce281cbcdaa5706f16f'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*_circularized_consensus.fasta"), emit: fasta
    //tuple val(meta), path("*_unpolished.log"),   emit: circlator_log
    path ("versions.yml"),                       emit: versions

    script:
    def container = task.container.toString() - "staphb/circlator@"
    """
    # If it finds an acceptable dnaA match, it rotates the circular sequence so that position 1 is at that gene.
    # If the gene is on the reverse strand, it can reverse-complement/reorient the contig so the selected gene is forward-facing.

    circlator fixstart $fasta ${meta.id}_circularized_consensus

    gzip ${meta.id}_circularized_consensus.fasta -c > ${meta.id}_circularized_consensus.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circlator: \$( circlator --version | sed -e "s/CIRCLATOR v//g" )
        circlator_container: ${container}
    END_VERSIONS
    """
}