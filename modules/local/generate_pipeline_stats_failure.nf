process GENERATE_PIPELINE_STATS_FAILURE {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.3.0 - MUST manually change below (line 56)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

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
    // terra=true sets paths for bc/wget for terra container paths
    def terra = params.terra ? "-2 terra" : ""
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "bash ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def srst_fullgenes_file = srst_fullgenes ? "-x $srst_fullgenes" : ""
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    # this runs with --mode CDC_PHEONIX or PHOENIX when SPAdes fails (creates contigs and not scaffolds)
    ${ica}pipeline_stats_writer.sh \\
        -a $raw_qc \\
        -b $fastp_total_qc \\
        -d ${prefix} \\
        -e $k2_trim_report \\
        -f $k2_trim_summary \\
        -g $krona_trimd \\
        -q $taxID \\
        $srst_fullgenes_file \\
        -5 $coverage \\
        $terra

    script_version=\$(${ica}pipeline_stats_writer.sh -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}
