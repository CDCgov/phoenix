process GENERATE_PIPELINE_STATS_FAILURE {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.0.0'

    input:
    tuple val(meta), path(raw_qc), \
    path(fastp_total_qc), \
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
    def prefix = task.ext.prefix ?: "${meta.id}"
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) {
        terra = ""
    } else if (params.terra==true) {
        terra = "-2 terra"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    pipeline_stats_writer.sh \\
        -a $raw_qc \\
        -b $fastp_total_qc \\
        -d ${prefix} \\
        -e $kraken2_trimd_report \\
        -f $kraken2_trimd_summary \\
        -g $krona_trimd \\
        -q $taxID \\
        $terra

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}