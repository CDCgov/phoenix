process BBDUK {
    tag "$meta.id"
    label 'process_medium'
    //v39.01
    container 'staphb/bbtools@sha256:161b0e1e198110b7edff8084ae9854d84eb32789d0fd62c7ced302078911c9d7'

    input:
    tuple val(meta), path(reads), val(fairy_outcome)
    path(contaminants)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_cleaned_1.fastq.gz out2=${prefix}_cleaned_2.fastq.gz"
    def contaminants_fa = contaminants ? "ref=$contaminants" : ''
    def maxmem = task.memory.toGiga()-(task.attempt*12) // keep heap mem low so and rest of mem is for java expansion.
    def container = task.container.toString() - "staphb/bbtools@"
    """
    maxmem=\$(echo \"$maxmem GB\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $args \\
        $contaminants_fa \\
        &> ${prefix}.bbduk.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
        bbmap_container: ${container}
    END_VERSIONS
    """
}