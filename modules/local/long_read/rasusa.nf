process RASUSA {
    tag "$meta"
    label 'process_high'
    container 'staphb/rasusa:0.7.0'
    //sha256:39642529ba1de30c1e52df8063c5381ffba793993a043df808dc2bdbaf655f78

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fastq.gz") , emit: fastq
    path "versions.yml"                       , emit: versions

    script:
    """
    rasusa --input $fastq -f 0.25 -o ${meta}_sub.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$( rasusa --version | sed -e "s/RASUSA v//g" )
    END_VERSIONS
    """

}
