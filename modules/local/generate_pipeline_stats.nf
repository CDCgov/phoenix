process GENERATE_PIPELINE_STATS {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 50)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

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
    path(gc_content),
    val(run_type)
    val(coverage)

    output:
    tuple val(meta), path('*.synopsis'), emit: pipeline_stats
    path("versions.yml")               , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def raw_qc_file            = raw_qc ? "--raw-read-counts $raw_qc" : "" // if raw_qc is null return "-a $raw_qc" else return ""
    def fastp_total_file       = fastp_total_qc ? "--total-read-counts $fastp_total_qc" : ""
    def k2_trim_report         = kraken2_trimd_report ? "--kraken2-trimd-report $kraken2_trimd_report" : ""
    def k2_trim_summary        = kraken2_trimd_summary ? "--kraken2-trimd-summary $kraken2_trimd_summary" : ""
    def krona_trim             = krona_trimd ? "--krona-trimd $krona_trimd" : ""
    
    def assembly_scaffolds_gz  = assembly_scaffolds ? "--assembly $assembly_scaffolds" : ""
    def filtered_assembly_gz   = filtered_assembly ? "--filtered-assembly $filtered_assembly" : ""
    def tax_file               = taxID ? "--taxid-file $taxID" : ""
    def quast_file             = quast_report ? "--quast-report $quast_report" : ""
    def mlst_file              = mlst ? "--mlst-file $mlst" : ""
    def assembly_ratio_file    = assembly_ratio ? "--assembly-ratio-file $assembly_ratio" : ""
    def gc_content_file        = gc_content ? "--gc-content-file $gc_content" : ""
    def fastANI_formatted_file = fastANI_formatted ? "--formatted-fastani $fastANI_formatted" : ""
    def k2_asmbld_report       = kraken2_asmbld_report ? "--kraken2-asmbld-report $kraken2_asmbld_report" : ""
    def k2_asmbld_summary      = kraken2_asmbld_summary ? "--kraken2-asmbled-summary $kraken2_asmbld_summary" : ""
    def krona_assembled        = krona_asmbld ? "--krona-asmbld $krona_asmbld" : ""
    def k2_wtasmbld_report     = kraken2_weighted_report ? "--kraken2-weighted-report $kraken2_weighted_report" : ""
    def k2_wtasmbld_summary    = kraken2_weighted_summary ? "--kraken2-weighted-summary $kraken2_weighted_summary" : ""
    def krona_wtasmbld         = krona_weighted ? "--krona-weighted $krona_weighted" : ""
    def ar_gamma_file          = ar_gamma ? "--gamma-ar $ar_gamma" : ""
    def amr_file               = amr_report ? "--amr-file $amr_report" : ""
    def pf_gamma_file          = pf_gamma ? "--gamma-replicon $pf_gamma" : ""
    def hv_gamma_file          = hv_gamma ? "--gamma-hv $hv_gamma" : ""
    def busco_summary          = busco_specific_short_summary ? "--busco-summary $busco_specific_short_summary" : ""
    def srst_fullgenes_file    = srst_fullgenes ? "--srst2-file $srst_fullgenes" : ""
    def extended_qc            = busco_specific_short_summary ? "--cdc-phoenix-mode" : ""
    def scaffold_flag          = (run_type == "SCAFFOLDS" || run_type == "CDC_SCAFFOLDS") ? "--assembly-only" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}pipeline_stats_writer.py \\
        $raw_qc_file \\
        $fastp_total_file \\
        --sample-name ${meta.id} \\
        $k2_trim_report \\
        $k2_trim_summary \\
        $krona_trim \\
        $assembly_scaffolds_gz \\
        $filtered_assembly_gz \\
        $k2_asmbld_report \\
        $k2_asmbld_summary \\
        $krona_assembled \\
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
        --coverage $coverage \\
        $extended_qc \\
        $scaffold_flag 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \$(${ica}pipeline_stats_writer.py -V)
    END_VERSIONS
    """
}
