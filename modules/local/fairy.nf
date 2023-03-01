process FAIRY {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_FAIry_synopsis.txt")    , optional:true,        emit: fairy_results

    script:
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Asses the FASTQ files for corruption
    echo "FASTQ_File:" ${prefix}_R1 >> ${prefix}_R1_FAIry.txt
    zcat ${reads[0]} | tail >> ${prefix}_R1_FAIry.txt

    echo "FASTQ_File:" ${prefix}_R2 >> ${prefix}_R2_FAIry.txt
    zcat ${reads[1]} | tail >> ${prefix}_R2_FAIry.txt
    
    if grep -Fxq "error" ${prefix}_R1_FAIry.txt 
        then echo ${prefix}_R1_FAIry.txt "is corrupted." >> ${prefix}_FAIry_synopsis.txt
        else echo ${prefix}_R1_FAIry.txt "is not corrupted." >> ${prefix}_FAIry_synopsis.txt
    fi

    if grep -Fxq "error" ${prefix}_R2_FAIry.txt 
        then echo ${prefix}_R2_FAIry.txt "is corrupted." >> ${prefix}_FAIry_synopsis.txt
        else echo ${prefix}_R2_FAIry.txt "is not corrupted." >> ${prefix}_FAIry_synopsis.txt
    fi
    """
}
