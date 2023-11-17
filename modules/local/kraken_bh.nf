process KRAKEN_BEST_HIT {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.0.0'

    input:
    tuple val(meta), path(kraken_summary), path(count_file) //[-q count_file (reads or congtigs)] so quast report for assembled or output of GATHERING_READ_QC_STATS for trimmed
    val(kraken_type) //weighted, trimmmed or assembled

    output:
    tuple val(meta), path('*.top_kraken_hit.txt'), emit:ksummary
    path("versions.yml")           , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) { terra = ""} 
    else if (params.terra==true) { terra = "-t terra" }
    else { error "Please set params.terra to either \"true\" or \"false\"" }
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables

    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}kraken2_best_hit.sh -i $kraken_summary -q $count_file -n ${prefix} $terra

    mv ${prefix}.summary.txt ${prefix}.kraken2_${kraken_type}.top_kraken_hit.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
        \$(${ica}kraken2_best_hit.sh -V)
    END_VERSIONS
    """
}