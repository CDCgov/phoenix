process MEDAKA {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/medaka:1.2.0'
    //sha256:3c84e7f69ada219a88ef06a10f187080fd5fab02d025f6eb048da2ad800186c2

    input:
    tuple val(meta), path(fasta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_medaka_consensus.fasta.gz"), emit: fasta_fin
    path ("versions.yml"),                                         emit: versions

    script:
    """
    medaka_consensus -i ${fastq} -d ${fasta} -o ./ -t $task.cpus

    cp consensus.fasta ${meta.id}_medaka_consensus.fasta

    gzip ${meta.id}_medaka_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version | sed -e "s/medaka//g" )
    END_VERSIONS
    """
}