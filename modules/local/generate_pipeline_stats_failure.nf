process GENERATE_PIPELINE_STATS_FAILURE {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 56)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta), path(raw_qc), \
    path(fastp_total_qc), \
    path(srst_fullgenes), \
    path(k2_trim_report), \
    path(krona_trimd), \
    path(k2_trim_summary), \
    path(taxID)
    val(coverage)

    output:
    tuple val(meta), path('*.synopsis'), emit: pipeline_stats
    path("versions.yml")               , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def srst_fullgenes_file = srst_fullgenes ? "--srst2-file $srst_fullgenes" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    # this runs with --mode CDC_PHEONIX or PHOENIX when SPAdes fails (creates contigs and not scaffolds)
    ${ica}pipeline_stats_writer.py \\
        --raw-read-counts $raw_qc \\
        --total-read-counts $fastp_total_qc \\
        --sample-name ${prefix} \\
        --kraken2-trimd-report $k2_trim_report \\
        --kraken2-trimd-summary $k2_trim_summary \\
        --krona-trimd $krona_trimd \\
        --taxid-file $taxID \\
        $srst_fullgenes_file \\
        --coverage $coverage

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \$(${ica}pipeline_stats_writer.py -V)
    END_VERSIONS
    """
}
