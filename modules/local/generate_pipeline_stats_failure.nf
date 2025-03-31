process GENERATE_PIPELINE_STATS_FAILURE {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 35)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

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
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) { terra = ""} 
    else if (params.terra==true) { terra = "-2 terra" }
    else { error "Please set params.terra to either \"true\" or \"false\"" }
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}pipeline_stats_writer.sh \\
        -a $raw_qc \\
        -b $fastp_total_qc \\
        -d ${prefix} \\
        -e $kraken2_trimd_report \\
        -f $kraken2_trimd_summary \\
        -g $krona_trimd \\
        -q $taxID \\
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