process GENERATE_PIPELINE_STATS_EXQC {
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
    path(assembly_scaffolds), \
    path(filtered_assembly), \
    path(mlst_file), \
    path(gamma_HV), \
    path(gamma_AR), \
    path(gamma_replicon), \
    path(quast_report), \
    path(busco_specific_short_summary), \
    path(kraken2_asmbld_report), \
    path(krona_asmbld), \
    path(kraken2_asmbled_summary), \
    path(krona_weighted), \
    path(kraken2_weighted_report), \
    path(kraken2_weighted_summary), \
    path(taxID), \
    path(fastANI_formatted_file), \
    path(assembly_ratio_file), \
    path(amr_file), \
    path(gc_content), \
    path(nanostat)
    val(coverage)

    output:
    tuple val(meta), path('*.synopsis'), emit: pipeline_stats
    path("versions.yml")               , emit: versions

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
    def raw             = raw_qc ? "-a $raw_qc" : "" // if raw_qc is null return "-a $raw_qc" else return ""
    def fastp_total     = fastp_total_qc ? "-b $fastp_total_qc" : ""
    def k2_trim_report  = kraken2_trimd_report ? "-e $kraken2_trimd_report" : ""
    def k2_trim_summary = kraken2_trimd_summary ? "-f $kraken2_trimd_summary" : ""
    def krona_trim      = krona_trimd ? "-g $krona_trimd" : ""
    def srst_file       = srst_fullgenes_file ? "-x $srst_fullgenes_file" : ""
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}pipeline_stats_writer.sh \\
        $raw \\
        $fastp_total \\
        -c $gc_content \\
        -d ${prefix} \\
        $k2_trim_report \\
        $k2_trim_summary \\
        $krona_trim \\
        -h $assembly_scaffolds \\
        -i $filtered_assembly \\
        -j $kraken2_asmbld_report \\
        -k $kraken2_asmbled_summary \\
        -l $krona_asmbld \\
        -m $kraken2_weighted_report \\
        -n $kraken2_weighted_summary \\
        -o $krona_weighted \\
        -p $quast_report \\
        -q $taxID \\
        -r $assembly_ratio_file \\
        -s $busco_specific_short_summary \\
        -t $fastANI_formatted_file \\
        -u $gamma_AR \\
        -v $gamma_replicon \\
        -w $gamma_HV \\
        $srst_file \\
        -y $mlst_file \\
        -4 $amr_file \\
        -5 $coverage \\
        -3 \\
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
