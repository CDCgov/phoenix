process FLYE {
    tag "$meta"
    label 'process_high'
    errorStrategy 'ignore'
    container 'staphb/flye:2.9.3'
    //sha256:bcca8c41798c6f4b3ce334a2faa9a010859fb2bbb2e8992a1b0abb84a8fa0c4c

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta}/*"), emit: dir
    tuple val(meta), path("${meta}/*.fasta"), emit: fasta
    tuple val(meta), path("${meta}/*.txt"), emit: assembly
    tuple val(meta), path("${meta}/*.gfa"), emit: gfa
    tuple val(meta), path("${meta}/*.tsv"), emit: flye_stat
    path "versions.yml"                       , emit: versions

    script:
    """
    flye --nano-hq $fastq -o ${meta} --iterations 1 -t 32
    awk -v contig="contig" 'BEGIN {OFS=FS="\t"} NR==1 {gsub(/^[^\t]+/, contig)} {print}' ${meta}/*.txt | cut -f 1-3 > ${meta}/${meta}_mqc.tsv 
    awk 'BEGIN {FS = "\t"; OFS = "," } { print } ' ${meta}/${meta}_mqc.tsv > ${meta}/${meta}_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version | sed -e "s/FLYE v//g" )
    END_VERSIONS
    """

}
