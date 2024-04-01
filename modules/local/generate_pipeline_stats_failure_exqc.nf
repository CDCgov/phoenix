process GENERATE_PIPELINE_STATS_FAILURE_EXQC {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 56)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

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
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script = params.ica ? "bash ${params.ica_path}/pipeline_stats_writer.sh" : "pipeline_stats_writer.sh"
    def terra = params.terra ? "-2 terra" : ""
    """
    # this runs with -entry CDC_PHEONIX when SPAdes fails (creates contigs and not scaffolds)
    ${script} \\
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

    script_version=\$(${script} -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}
