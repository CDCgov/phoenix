process RENAME_SRA_FASTA {
    tag "${meta.id}"
    label 'process_low'
    // base_v2.3.0 - MUST manually change below (line 16)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(reads)

    output:
    path("*_R*_001.fastq.gz"), emit: renamed_reads // we don't want the SRR.fastq just the forward and reverse
    path("versions.yml"),      emit: versions

    script:
    def srr_num = reads[0].toString() - "_1.fastq.gz" // this is the SRR number
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    mv ${srr_num}_1.fastq.gz ${meta.id}_R1_001.fastq.gz
    mv ${srr_num}_2.fastq.gz ${meta.id}_R2_001.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}