process POLYPOLISH {
    tag "${meta}"
    label 'process_high'
    container 'quay.io/staphb/polypolish'
    errorStrategy 'ignore'
    

    input:
    tuple val(meta), path (fasta), file(sam)

    output:
    tuple val(meta), path("${meta.id}_polished_consensus.fasta.gz"), emit: assembly
    path "versions.yml",                                             emit: versions

    script:
    """
    polypolish filter --in1 ${meta.id}_1.sam --in2 ${meta.id}_2.sam --out1 ${meta.id}_filtered_1.sam --out2 ${meta.id}_filtered_2.sam
    polypolish polish --careful $fasta ${meta.id}_filtered_1.sam ${meta.id}_filtered_2.sam > ${meta.id}_polished_consensus.fasta
    #header.sh ${meta.id}_polished_consensus.fasta

    #gzip file for down stream process
    gzip --force ${meta.id}_polished_consensus.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        polypolish: \$( polypolish --version | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}
