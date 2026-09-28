process RASUSA {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/rasusa@sha256:c16b8154a90e9dbe4037d5dc1af5332e3ab8d18ebe2753fafc02964fe3daed27'

    input:
    tuple val(meta), path(reads), path(estimation)
    val(depth)

    output:
    tuple val(meta), path("*_${depth}X.fastq.gz"), emit: subfastq
    path ("versions.yml"),                         emit: versions

    script:
    def container = task.container.toString() - "staphb/rasusa@"
    """
    echo "FASTQ: $reads" > error.txt
    GENOME_SIZE=\$(<$estimation)
    echo "Estimated size: \$GENOME_SIZE" >> error.txt

    rasusa reads --genome-size \$GENOME_SIZE --coverage $depth -o ${meta.id}_${depth}X.fastq.gz --seed 42 $reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(rasusa --version | sed -e "s/rasusa //g")
        rasusa_container: ${container}
    END_VERSIONS
    """

}
