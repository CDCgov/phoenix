process FASTANI {
    tag "$meta.id"
    label 'process_medium'
    beforeScript 'ulimit -s unlimited'
    container 'staphb/fastani:1.33'
    /*container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.33--h0fdf51a_0' :
        'quay.io/biocontainers/fastani:1.33--h0fdf51a_0' }"*/

    input:
    tuple val(meta), path(query), path(reference), path(reference_files)
    

    output:
    tuple val(meta), path("*.ani.txt"), emit: ani
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastANI \\
        -q $query \\
        --rl $reference \\
        -o ${prefix}.ani.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}