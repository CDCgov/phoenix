process POLYPOLISH {
    tag "${meta}"
    label 'process_high'
    container 'staphb/polypolish@sha256:46557d6b6c0373a95ec34c0b90c8a10a25ae6b64e601f073474494a4fedc6b10'

    input:
    tuple val(meta), path (fasta), file(sam)

    output:
    tuple val(meta), path("${meta.id}_polished_consensus.fasta.gz"), emit: assembly  // fasta.gz for emits in hybrid.nf and modules.config
    tuple val(meta), path("${meta.id}_polished_consensus.fasta"),    emit: assembly_fasta  // fasta for input of plasmid_characterization subworkflow 
    path "versions.yml",                                             emit: versions

    script:
    def container = task.container.toString() - "staphb/polypolish@"
    """
    polypolish filter --in1 ${meta.id}_1.sam --in2 ${meta.id}_2.sam --out1 ${meta.id}_filtered_1.sam --out2 ${meta.id}_filtered_2.sam
    polypolish polish --careful $fasta ${meta.id}_filtered_1.sam ${meta.id}_filtered_2.sam > ${meta.id}_polished_consensus.fasta
    #header.sh ${meta.id}_polished_consensus.fasta

    #gzip file for down stream process
    gzip --force -k ${meta.id}_polished_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$( polypolish --version | cut -f 2 -d ' ' )
        polypolish_container: ${container}
    END_VERSIONS
    """
}
