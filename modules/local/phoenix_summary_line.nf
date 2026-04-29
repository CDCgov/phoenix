process CREATE_SUMMARY_LINE {
    tag "${meta.id}"   // <-- closure for tag
    label 'process_single'
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

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
        path(shigapass),
        val(old_software_version)
    val(extended_qc)
    val(phx_version)

    output:
    tuple val(meta), path('*_summaryline.tsv'), emit: line_summary
    path("versions.yml"),                       emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    
    // We use .peek() or checking name/size to ensure these only create flags if files actually exist
    def extended_qc_arg        = extended_qc ? "--extended_qc" : ""
    
    // Using Groovy's ability to check if a file path is actually a valid file name
    def trim_ksummary_file     = (trimd_ksummary && !(trimd_ksummary.toString() == "[]")) ? "-k $trimd_ksummary" : ""
    def wtasmbld_ksummary_file = (wtasmbld_ksummary && !(wtasmbld_ksummary.toString() == "[]")) ? "--kraken_wtasmbld ${wtasmbld_ksummary}" : ""
    def fastani_file           = (fastani && !(fastani.toString() == "[]")) ? "-f $fastani" : ""
    def quast_file             = (quast_report && !(quast_report.toString() == "[]")) ? "-q $quast_report" : ""
    def ar_gamma_file          = (ar_gamma && !(ar_gamma.toString() == "[]")) ? "-a $ar_gamma" : ""
    def amr_file               = (amr_report && !(amr_report.toString() == "[]")) ? "-u $amr_report" : ""
    def pf_gamma_file          = (pf_gamma && !(pf_gamma.toString() == "[]")) ? "-p $pf_gamma" : ""
    def hv_gamma_file          = (hv_gamma && !(hv_gamma.toString() == "[]")) ? "-v $hv_gamma" : ""
    def ratio_file             = (assembly_ratio && !(assembly_ratio.toString() == "[]")) ? "-r $assembly_ratio" : ""
    def mlst_file              = (mlst && !(mlst.toString() == "[]")) ? "-m $mlst" : ""
    def fastp_file             = (fastp_total_qc && !(fastp_total_qc.toString() == "[]")) ? "-t $fastp_total_qc" : ""
    // Check if we actually got a version string
    def old_software_arg = (old_software_version && old_software_version != "Unknown") ? "--old_phoenix_version ${old_software_version}" : ""
    
    // Shigapass needs careful joining because it might be a list or a single path
    def shigapass_file = (shigapass && shigapass.toString() != "[]") ? "--shigapass ${shigapass instanceof List ? shigapass.join(' ') : shigapass}" : ""

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
        $old_software_arg \\
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
