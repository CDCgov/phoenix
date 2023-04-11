process FASTP_SINGLES {
    tag "$meta.id"
    label 'process_low'
    container 'staphb/fastp:0.23.2'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.singles.fastq.gz')  , emit: reads
    tuple val(meta), path('*.json')              , emit: json
    tuple val(meta), path('*.html')              , emit: html
    tuple val(meta), path('*.log')               , emit: log
    path "versions.yml"                          , emit: versions
    tuple val(meta), path('*.merged.fastq.gz')   , optional:true, emit: reads_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ ! -s ${reads[0]} ]] && [[ ! -s ${reads[1]} ]]; then
        # Both are empty, do nothing???
        :
    elif [[ ! -s ${reads[0]} ]]; then
        # Only R1 is empty, run on R2 only
        zcat ${reads[1]} > ${prefix}.cat_singles.fastq
        gzip ${prefix}.cat_singles.fastq
    elif [[ ! -s ${reads[1]} ]]; then
        # Only R2 is empty, run on R1 only
        zcat ${reads[0]} > ${prefix}.cat_singles.fastq
        gzip ${prefix}.cat_singles.fastq
    else
        # Both reads have contents
        zcat ${reads[0]} ${reads[1]} > ${prefix}.cat_singles.fastq
        gzip ${prefix}.cat_singles.fastq
    fi
    # Possibly will need to catch when in1 is empty, dont know how fastp handles it, but right now we need to many of its standard outputs
        fastp \\
            --in1 ${prefix}.cat_singles.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}_singles.fastp.json \\
            --html ${prefix}_singles.fastp.html \\
            --out1 ${prefix}.singles.fastq.gz \\
            $args \\
            2> ${prefix}.fastp.log


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}