process GENERATE_PIPELINE_STATS_FAILURE_EXQC {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(raw_qc), \
    path(fastp_total_qc), \
    path(srst_fullgenes_file), \
    path(kraken2_trimd_report), \
    path(krona_trimd), \
    path(kraken2_trimd_summary), \
    path(taxID), \
    val(spades_outcome)
    val(coverage)

    output:
    tuple val(meta), path('*.synopsis'), emit: pipeline_stats
    path("versions.yml")               , emit: versions

    when:
    "${spades_outcome[0]}" == "run_failure" || "${spades_outcome[1]}" == "no_scaffolds" || "${spades_outcome[2]}" == "no_contigs"

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) { terra = ""} 
    else if (params.terra==true) { terra = "-2 terra" }
    else { error "Please set params.terra to either \"true\" or \"false\"" }
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    # this runs with -entry CDC_PHEONIX when SPAdes fails (creates contigs and not scaffolds)
    ${ica}pipeline_stats_writer.sh \\
        -a $raw_qc \\
        -b $fastp_total_qc \\
        -d ${prefix} \\
        -e $kraken2_trimd_report \\
        -f $kraken2_trimd_summary \\
        -g $krona_trimd \\
        -q $taxID \\
        -x $srst_fullgenes_file \\
        -5 $coverage \\
        $terra

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
