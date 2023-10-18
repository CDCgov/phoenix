process QUAST {
    tag "$meta.id"
    label 'process_low'
    container 'staphb/quast:5.0.2'

    input:
    tuple val(meta), path(consensus), val(fairy_outcome)

    output:
    tuple val(meta), path("quast")        , emit: results
    tuple val(meta), path('*.tsv')        , emit: report_tsv
    path "versions.yml"                   , emit: versions

    when:
    //if the files are not corrupt and there are equal number of reads in each file then run bbduk
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    quast.py \\
        --output-dir quast \\
        --threads $task.cpus \\
        $args \\
        $consensus

    mv quast/report.tsv ./${prefix}_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | grep "QUAST" | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
