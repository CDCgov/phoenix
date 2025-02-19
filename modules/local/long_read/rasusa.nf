process RASUSA {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/rasusa:0.7.0'
    //sha256:39642529ba1de30c1e52df8063c5381ffba793993a043df808dc2bdbaf655f78

    input:
    tuple val(meta), path(trimmed_fastq)

    output:
    tuple val(meta), path("*_sub.fastq.gz"), emit: fastq
    path ("versions.yml"),                   emit: versions

    script:
    """
    rasusa --input ${trimmed_fastq} --frac 0.25 -o ${meta.id}_sub.fastq.gz --seed 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$( rasusa --version | sed -e "s/RASUSA v//g" )
    END_VERSIONS
    """

}
