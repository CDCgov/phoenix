process GENERATE_PIPELINE_STATS {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 50)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f7cb3aa4e3324cab43d8635be17da8ae15f62e39d380acda844d1c9deef69c60'

    input:
    tuple val(meta), path(raw_qc), \
    path(fastp_total_qc), \
    path(srst_fullgenes), \
    path(kraken2_trimd_report), \
    path(krona_trimd), \
    path(kraken2_trimd_summary), \
    path(assembly_scaffolds), \
    path(filtered_assembly), \
    path(mlst), \
    path(hv_gamma), \
    path(ar_gamma), \
    path(pf_gamma), \
    path(quast_report), \
    path(busco_specific_short_summary), \
    path(kraken2_asmbld_report), \
    path(krona_asmbld), \
    path(kraken2_asmbld_summary), \
    path(krona_weighted), \
    path(kraken2_weighted_report), \
    path(kraken2_weighted_summary), \
    path(taxID), \
    path(fastANI_formatted), \
    path(assembly_ratio), \
    path(amr_report), \
    path(gc_content)
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
    def raw_qc_file            = raw_qc ? "-a $raw_qc" : "" // if raw_qc is null return "-a $raw_qc" else return ""
    def fastp_total_file       = fastp_total_qc ? "-b $fastp_total_qc" : ""
    def k2_trim_report         = kraken2_trimd_report ? "-e $kraken2_trimd_report" : ""
    def k2_trim_summary        = kraken2_trimd_summary ? "-f $kraken2_trimd_summary" : ""
    def krona_trim             = krona_trimd ? "-g $krona_trimd" : ""
    
    def assembly_scaffolds_gz  = krona_trimd ? "-h $assembly_scaffolds" : ""
    def filtered_assembly_gz   = krona_trimd ? "-i $filtered_assembly" : ""
    def tax_file               = taxID ? "-q $taxID" : ""
    def quast_file             = quast_report ? "-p $quast_report" : ""
    def mlst_file              = mlst ? "-y $mlst" : ""
    def assembly_ratio_file    = assembly_ratio ? "-r $assembly_ratio" : ""
    def gc_content_file        = gc_content ? "-c $gc_content" : ""
    def fastANI_formatted_file = fastANI_formatted ? "-t $fastANI_formatted" : ""

    def k2_asmbld_report       = kraken2_asmbld_report ? "-j $kraken2_asmbld_report" : ""
    def k2_asmbld_summary      = kraken2_asmbld_summary ? "-k $kraken2_asmbld_summary" : ""
    def krona_asmbld           = krona_asmbld ? "-l $krona_asmbld" : ""
    def k2_wtasmbld_report     = kraken2_weighted_report ? "-m $kraken2_weighted_report" : ""
    def k2_wtasmbld_summary    = kraken2_weighted_summary ? "-n $kraken2_weighted_summary" : ""
    def krona_wtasmbld         = krona_weighted ? "-o $krona_weighted" : ""
    
    def ar_gamma_file          = ar_gamma ? "-u $ar_gamma" : ""
    def amr_file               = amr_report ? "-4 $amr_report" : ""
    def pf_gamma_file          = ar_gamma ? "-v $pf_gamma" : ""
    def hv_gamma_file          = hv_gamma ? "-w $hv_gamma" : ""
    def busco_summary          = busco_specific_short_summary ? "-s $busco_specific_short_summary" : ""
    def srst_fullgenes_file    = srst_fullgenes ? "-x $srst_fullgenes" : ""
    def extended_qc            = busco_specific_short_summary ? "-3" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}pipeline_stats_writer.sh \\
        $raw_qc_file \\
        $fastp_total_file \\
        -d ${meta.id} \\
        $k2_trim_report \\
        $k2_trim_summary \\
        $krona_trim \\
        $assembly_scaffolds_gz \\
        $filtered_assembly_gz \\
        $k2_asmbld_report \\
        $k2_asmbld_summary \\
        $krona_asmbld \\
        $k2_wtasmbld_report \\
        $k2_wtasmbld_summary \\
        $krona_wtasmbld \\
        $quast_file \\
        $tax_file \\
        $gc_content_file \\
        $assembly_ratio_file \\
        $busco_summary \\
        $fastANI_formatted_file \\
        $ar_gamma_file \\
        $pf_gamma_file \\
        $hv_gamma_file \\
        $srst_fullgenes_file \\
        $mlst_file \\
        $amr_file \\
        -5 $coverage \\
        $extended_qc \\
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
