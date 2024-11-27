process MEDAKA {
    tag "$meta"
    label 'process_high'
    errorStrategy 'ignore'
    container 'staphb/medaka:1.2.0'
    //sha256:3c84e7f69ada219a88ef06a10f187080fd5fab02d025f6eb048da2ad800186c2

    input:
    tuple val(meta), path(fasta), path(fastq)

    output:
    path "${meta}/",                                                 emit: directory
    tuple val(meta), path("${meta}/${meta}_medaka_consensus.fasta"), emit: fasta_fin
    path "versions.yml"                       , emit: versions

script:
  """

    medaka_consensus -i ${fastq} -d ${fasta} -o ${meta} -t 32
    cp ${meta}/consensus.fasta ${meta}/${meta}_medaka_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version | sed -e "s/MEDAKA v//g" )
    END_VERSIONS
  """
}