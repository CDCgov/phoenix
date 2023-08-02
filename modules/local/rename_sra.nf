process RENAME_SRA_FASTA {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.0.2'

    input:
    tuple val(meta), path(reads)

    output:
    path("*_R*_001.fastq.gz"), emit: renamed_reads // we don't want the SRR.fastq just the forward and reverse
    path("versions.yml"),      emit: versions

    script:
    def srr_num = reads[0].toString() - "_1.fastq.gz" // this is the SRR number
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    mv ${srr_num}_1.fastq.gz ${meta.id}_R1_001.fastq.gz
    mv ${srr_num}_2.fastq.gz ${meta.id}_R2_001.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}