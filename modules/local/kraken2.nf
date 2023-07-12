process KRAKEN2_KRAKEN2 {
    tag "$meta.id"
    label 'process_high'
    container 'staphb/kraken2:2.1.2-no-db'

    input:
    tuple val(meta), path(reads), path(db)
    val(kraken_type) //weighted, trimmmed or assembled
    val(save_output_fastqs)
    val(save_reads_assignment)

    output:
    tuple val(meta), path('*classified*')     , optional:true, emit: classified_reads_fastq
    tuple val(meta), path('*unclassified*')   , optional:true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads*'), optional:true, emit: classified_reads_assignment
    tuple val(meta), path('*report.txt')                     , emit: report
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                       = task.ext.args ?: ''
    def prefix                     = task.ext.prefix ?: "${meta.id}"
    def paired                     = meta.single_end ? "" : "--paired"
    def classified                 = meta.single_end ? "${prefix}.classified.fasta"   : "${prefix}.classified#.fasta"
    def unclassified               = meta.single_end ? "${prefix}.unclassified.fasta" : "${prefix}.unclassified#.fasta"
    def classified_command         = save_output_fastqs ? "--classified-out ${classified}" : ""
    def unclassified_command       = save_output_fastqs ? "--unclassified-out ${unclassified}" : ""
    def readclassification_command = save_reads_assignment ? "--output ${prefix}.kraken2_${kraken_type}.classifiedreads.txt" : ""
    def compress_reads_command     = save_output_fastqs ? "gzip *.fasta" : ""

    """
    kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${prefix}.kraken2_${kraken_type}.report.txt \\
        --gzip-compressed \\
        $unclassified_command \\
        $classified_command \\
        $readclassification_command \\
        $paired \\
        $args \\
        $reads

    $compress_reads_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
        kraken2db: \$(echo \$(basename "${params.kraken2db}"))
    END_VERSIONS
    """
}
