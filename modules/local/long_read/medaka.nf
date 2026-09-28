process MEDAKA {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/medaka@sha256:be8deb0b6e901a0665e357ff3e7c16f108a4df96da2813ca44ff2fcd36a11b48'

    input:
    tuple val(meta), path(fasta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}_consensus.fasta"),    emit: fasta //for plasmid_characterization subworkflow
    tuple val(meta), path("${meta.id}_consensus.fasta.gz"), emit: fasta_gz  // for SCAFFOLDS_EXTERNAL workflow in main.nf
    path ("versions.yml"),                                  emit: versions

    script:
    def container = task.container.toString() - "staphb/medaka@"
    """
    medaka_consensus -i ${fastq} -d ${fasta} -o ${meta.id} -t 16 --bacteria

    cp ${meta.id}/consensus.fasta ${meta.id}_consensus.fasta

    gzip ${meta.id}_consensus.fasta -c > ${meta.id}_consensus.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version | sed -e "s/medaka//g" )
        medaka_container: ${container}
    END_VERSIONS
    """
}