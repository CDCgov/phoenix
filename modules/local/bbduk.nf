process BBDUK {
    tag "$meta.id"
    label 'process_medium'
    container 'staphb/bbtools:39.01' 

    input:
    tuple val(meta), path(reads)
    path contaminants 
    tuple val(meta), path(outcome)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_cleaned_1.fastq.gz out2=${prefix}_cleaned_2.fastq.gz"
    def contaminants_fa = contaminants ? "ref=$contaminants" : ''
    def maxmem = task.memory.toGiga()-(task.attempt*12) // keep heap mem low so and rest of mem is for java expansion.
    def checker = "{outcome}".minus("_raw_read_counts.txt")
    """
    if grep -Fxq "PASS" ${outcome} && [[ $raw = $checker* ]]
    then
        maxmem=\$(echo \"$maxmem GB\"| sed 's/ GB/g/g')
        bbduk.sh \\
            -Xmx\$maxmem \\
            $raw \\
            $trimmed \\
            threads=$task.cpus \\
            $args \\
            $contaminants_fa \\
            &> ${prefix}.bbduk.log
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
    END_VERSIONS
    """
}
