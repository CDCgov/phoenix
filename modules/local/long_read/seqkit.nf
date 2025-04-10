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

    

    script:
    """
    seqkit stats -Ta ${reads} > ${meta.id}_raw.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit --version | sed -e "s/seqkit //g")
    END_VERSIONS
    
    """
}
