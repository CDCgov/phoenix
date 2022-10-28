process BBMAP_REFORMAT {
    tag "$meta.id"
    label 'process_medium'
    container 'staphb/bbtools:38.96'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*filtered.scaffolds.fa.gz')   , emit: reads
    tuple val(meta), path('*.log')                       , emit: log
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = "in=${reads[0]}"
    def trimmed  = "out=${prefix}.filtered.scaffolds.fa.gz"
    def minlength = params.minlength
    def maxmem = task.memory.toGiga()-8 // keep heap mem low so and rest of space to java expansion.
    """
    maxmem=\$(echo \"$maxmem GB\"| sed 's/ GB/g/g')
    reformat.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $args \\
        minlength=$minlength \\
        &> ${prefix}.bbmap_filtered.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
