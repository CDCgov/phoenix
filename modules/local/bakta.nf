process BAKTA_BAKTA {
    tag "$meta.id"
    label 'process_medium'
    //
    container 'staphb/bakta@sha256:26710927a05791141dcfedf58b58c1f11cc919b199315dbd85ed07654342c694'

    input:
    tuple val(meta), path(fasta)
    path db
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("${prefix}.embl")             , emit: embl
    tuple val(meta), path("${prefix}.faa")              , emit: faa
    tuple val(meta), path("${prefix}.ffn")              , emit: ffn
    tuple val(meta), path("${prefix}.fna")              , emit: fna
    tuple val(meta), path("${prefix}.gbff")             , emit: gbff
    tuple val(meta), path("${prefix}.gff3")             , emit: gff
    tuple val(meta), path("${prefix}.hypotheticals.tsv"), emit: hypotheticals_tsv
    tuple val(meta), path("${prefix}.hypotheticals.faa"), emit: hypotheticals_faa
    tuple val(meta), path("${prefix}.tsv")              , emit: tsv
    tuple val(meta), path("${prefix}.txt")              , emit: txt
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_tf = prodigal_tf ? "--prodigal-tf ${prodigal_tf[0]}" : ""
    """
    bakta \\
        $fasta \\
        $args \\
        --threads $task.cpus \\
        --prefix $prefix \\
        $proteins_opt \\
        $prodigal_tf \\
        --db $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(echo \$(bakta --version) 2>&1 | cut -f '2' -d ' ')
        amrfinderplus_db_version: \$(head $db/amrfinderplus-db/latest/version.txt)
    END_VERSIONS
    """
}