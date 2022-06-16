process KRAKEN_BEST_HIT {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(kraken_report), path(count_file) //[-q count_file (reads or congtigs)] so quast report for assembled or output of GATHERING_READ_QC_STATS for trimmed
    val(kraken_type) //weighted, trimmmed or assembled

    output:
    tuple val(meta), path('*_summary.txt'), emit:ksummary

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kraken2_best_hit.sh -i $kraken_report -q $count_file -n ${prefix}

    mv ${prefix}.summary.txt ${prefix}.${kraken_type}_summary.txt
    """
}