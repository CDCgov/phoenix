process RENAME_SRA_FASTA {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(reads)

    output:
    path("*_*.fastq.gz")   , emit: renamed_reads // we don't want the SRR.fastq just the forward and reverse

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """ 
    mv ${prefix}_1.fastq.gz ${prefix}_R1_001.fastq.gz
    mv ${prefix}_2.fastq.gz ${prefix}_R2_001.fastq.gz
    """
}