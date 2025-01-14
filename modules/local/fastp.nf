process FASTP {
    tag "$meta.id"
    label 'process_medium'
    // v0.24.0
    container 'staphb/fastp@sha256:22ef1eedf8a5529e9ebab2b4aa77308d7f4235a0f11a4f7190e28d683b156702'

    input:
    tuple val(meta), path(reads)
    val(save_trimmed_fail)
    val(save_merged)

    output:
    tuple val(meta), path('*.trim.fastq.gz')  , optional:true, emit: reads
    tuple val(meta), path('*.json')           , emit: json
    tuple val(meta), path('*.html')           , emit: html
    tuple val(meta), path('*.log')            , emit: log
    path "versions.yml"                       , emit: versions
    tuple val(meta), path('*.fail.fastq.gz')  , optional:true, emit: reads_fail
    tuple val(meta), path('*.merged.fastq.gz'), optional:true, emit: reads_merged

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def container = task.container.toString() - "staphb/fastp@"
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        def fail_fastq = save_trimmed_fail ? "--failed_out ${prefix}.fail.fastq.gz" : ''
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        fastp \\
            --in1 ${prefix}.fastq.gz \\
            --out1 ${prefix}.trim.fastq.gz \\
            --thread $task.cpus \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $args \\
            2> ${prefix}.fastp.log
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        def fail_fastq  = save_trimmed_fail ? "--unpaired1 ${prefix}_1.fail.fastq.gz --unpaired2 ${prefix}_2.fail.fastq.gz" : ''
        def merge_fastq = save_merged ? "-m --merged_out ${prefix}.merged.fastq.gz" : ''
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        fastp \\
            --in1 ${prefix}_1.fastq.gz \\
            --in2 ${prefix}_2.fastq.gz \\
            --out1 ${prefix}_1.trim.fastq.gz \\
            --out2 ${prefix}_2.trim.fastq.gz \\
            --json ${prefix}.fastp.json \\
            --html ${prefix}.fastp.html \\
            $fail_fastq \\
            $merge_fastq \\
            --thread $task.cpus \\
            --detect_adapter_for_pe \\
            $args \\
            2> ${prefix}.fastp.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
            fastp_container: ${container}
        END_VERSIONS
        """
    }
}