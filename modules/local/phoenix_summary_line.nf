process CREATE_SUMMARY_LINE {
    tag "${meta.id}"   // <-- closure for tag
    label 'process_single'
    container 'quay.io/jvhagey/phoenix@sha256:f7cb3aa4e3324cab43d8635be17da8ae15f62e39d380acda844d1c9deef69c60'

    input:
    tuple val(meta),
        path(fastp_total_qc),
        path(mlst),
        path(hv_gamma),
        path(ar_gamma),
        path(pf_gamma),
        path(quast_report),
        path(assembly_ratio),
        path(synopsis),
        path(taxonomy_file),
        path(trimd_ksummary),
        path(wtasmbld_ksummary),
        path(amr_report),
        path(fastani),
        path(shigapass)
    val(extended_qc)
    val(phx_version)

    output:
    tuple val(meta), path('*_summaryline.tsv'), emit: line_summary
    path("versions.yml"),                       emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    // allowing for some optional parameters for -entry SCAFFOLDS/CDC_SCAFFOLDS nothing should be passed.
    def extended_qc_arg        = extended_qc ? "--extended_qc" : ""  // only for spades failures
    def trim_ksummary_file     = trimd_ksummary ? "-k $trimd_ksummary" : ""
    def wtasmbld_ksummary_file = wtasmbld_ksummary ? "--kraken_wtasmbld ${wtasmbld_ksummary}" : ""
    def fastani_file           = fastani ? "-f $fastani" : ""
    def quast_file             = quast_report ? "-q $quast_report" : ""
    def ar_gamma_file          = ar_gamma ? "-a $ar_gamma" : ""
    def amr_file               = amr_report ? "-u $amr_report" : ""
    def pf_gamma_file          = pf_gamma ? "-p $pf_gamma" : ""
    def hv_gamma_file          = hv_gamma ? "-v $hv_gamma" : ""
    def ratio_file             = assembly_ratio ? "-r $assembly_ratio" : ""
    def mlst_file              = mlst ? "-m $mlst" : ""
    def fastp_file             = fastp_total_qc ? "-t $fastp_total_qc" : ""
    def shigapass_file = (shigapass && shigapass.size() > 0) ? "--shigapass ${shigapass.join(' ')}" : ''
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Phoenix_summary_line.py \\
        $quast_file \\
        $fastp_file \\
        $ar_gamma_file \\
        $hv_gamma_file \\
        $pf_gamma_file \\
        $ratio_file \\
        $mlst_file \\
        $amr_file \\
        -n ${meta.id} \\
        -s $synopsis \\
        -x $taxonomy_file \\
        $fastani_file \\
        $wtasmbld_ksummary_file \\
        $trim_ksummary_file \\
        $shigapass_file \\
        --phx_version $phx_version \\
        -o ${meta.id}_summaryline.tsv \\
        $extended_qc_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        Phoenix_summary_line.py: \$(${ica}Phoenix_summary_line.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
