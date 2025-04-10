process FLYE {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/flye:2.9.3'
    //sha256:bcca8c41798c6f4b3ce334a2faa9a010859fb2bbb2e8992a1b0abb84a8fa0c4c
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.txt"),   emit: assembly_info
    tuple val(meta), path("*.gfa"),   emit: gfa
    tuple val(meta), path("*.tsv"),   emit: flye_stat
    path ("versions.yml"),            emit: versions

    script:
    """
    flye --nano-hq $fastq -o . --iterations 1 --threads $task.cpus

    # rename output files to include the meta.id
    mv assembly_info.txt ${meta.id}.assembly_info.txt
    mv flye.log ${meta.id}.flye.log

    # Comment what this is for - pulling information out of assembly_info.txt, why not just leave it there?
    awk -v contig="contig" 'BEGIN {OFS=FS="\t"} NR==1 {gsub(/^[^\t]+/, contig)} {print}' *.txt | cut -f 1-3 > ${meta.id}_mqc.tsv 
    awk 'BEGIN {FS = "\t"; OFS = "," } { print } ' ${meta.id}_mqc.tsv > ${meta.id}_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version | sed -e "s/FLYE v//g" )
    END_VERSIONS
    """
}
