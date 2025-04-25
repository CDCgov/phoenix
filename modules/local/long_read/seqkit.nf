process RAWSTATS {
    tag "${meta.id}"
    label 'process_medium'
    errorStrategy 'ignore'
    container 'staphb/seqkit'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_raw.txt"),         emit: rawstats
    path ("versions.yml"),                      emit: versions
    tuple val(meta), path("${meta.id}_lr.fastq.gz"), emit: fastq_lr

    

    script:
    """
    seqkit stats -Ta ${reads} > ${meta.id}_raw.txt
    mv $reads ${meta.id}_lr.fastq.gz
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit --version | sed -e "s/seqkit //g")
    END_VERSIONS
    
    """
}
