process LRGE {
    tag "${meta.id}"
    label 'process_medium'
    //errorStrategy 'ignore'
    container 'quay.io/staphb/lrge'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(fastq_lr)

    output:
    tuple val(meta), path("*_size.txt") , emit: estimation
    path ("versions.yml"),                      emit: versions

    

    script:
    """
    lrge -t $task.cpus ${fastq_lr} > ${meta.id}_size.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LRGE: \$(LRGE --version | sed -e "s/LRGE//g")

    END_VERSIONS
    """
}
