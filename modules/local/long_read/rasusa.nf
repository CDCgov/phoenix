process RASUSA {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/rasusa:0.7.0'
    //sha256:39642529ba1de30c1e52df8063c5381ffba793993a043df808dc2bdbaf655f78

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(estimation)
    val depth

    output:
    tuple val(meta), path("*_${depth}X.fastq.gz"), emit: subfastq
    path ("versions.yml"),                   emit: versions

    script:
    """
    cat $estimation > error.txt
    echo "FASTQ: $reads" >> error.txt
    GENOME_SIZE=\$(<$estimation)
    echo "Estimated size: \$GENOME_SIZE" >> error.txt
    GENOME_SIZE=\$(<$estimation)
    rasusa --input ${reads} --genome-size \$GENOME_SIZE --coverage $depth -o ${meta.id}_${depth}X.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(rasusa --version | sed -e "s/rasusa //g")
  
    END_VERSIONS
    """

}
