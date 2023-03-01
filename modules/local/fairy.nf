process FAIRY {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_R1_FAIry.txt")    , optional:true,        emit: fairy_results_1
    tuple val(meta), path("*_R2_FAIry.txt")    , optional:true,        emit: fairy_results_2

    script:
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "FASTQ_File:" ${prefix}_R1 >> ${prefix}_R1_FAIry.txt
    zcat ${reads[0]} | tail >> ${prefix}_R1_FAIry.txt

    echo "FASTQ_File:" ${prefix}_R2 >> ${prefix}_R2_FAIry.txt
    zcat ${reads[1]} | tail >> ${prefix}_R2_FAIry.txt
    """
}
