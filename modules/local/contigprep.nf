process CONTIG_PREP {
    tag "$meta.id"
    label 'process_medium'
    container 'staphb/bbtools:39.01'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.contigs.fa.gz')   , emit: reads
    path('*.contigs.fa.gz')                    , emit: contigs
    tuple val(meta), path('*.log')             , emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.contigs.fa.gz" : "out1=${prefix}_2.contigs.fa.gz out2=${prefix}_4.contigs.fa.gz"
    def minlength = params.minlength
    """
    reformat.sh \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $args \\
        minlength=$minlength \\
        &> ${prefix}.reformat.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
