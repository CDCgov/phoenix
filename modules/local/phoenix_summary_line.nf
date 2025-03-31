process CREATE_SUMMARY_LINE {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    tuple val(meta), path(trimmed_qc_data_file), \
    path(mlst_file), \
    path(hypervirulence_gamma_file), \
    path(ar_gamma_file), \
    path(pf_gamma_file), \
    path(quast_report), \
    path(ratio_file), \
    path(synopsis), \
    path(taxonomy_file), \
    path(trimd_ksummary), \
    path(amr_report), \
    path(fastani)

    output:
    tuple val(meta), path('*_summaryline.tsv'), emit: line_summary
    path("versions.yml"),                       emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    // allowing for some optional parameters for -entry SCAFFOLDS/CDC_SCAFFOLDS nothing should be passed.
    def trimmed_qc_data = trimmed_qc_data_file ? "-t $trimmed_qc_data_file" : ""
    def trim_ksummary   = trimd_ksummary ? "-k $trimd_ksummary" : ""
    def fastani_file    = fastani ? "-f $fastani" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Phoenix_summary_line.py \\
        -q $quast_report \\
        $trimmed_qc_data \\
        -a $ar_gamma_file \\
        -v $hypervirulence_gamma_file \\
        -p $pf_gamma_file \\
        -r $ratio_file \\
        -m $mlst_file \\
        -u $amr_report \\
        -n ${prefix} \\
        -s $synopsis \\
        -x $taxonomy_file \\
        $fastani_file \\
        $trim_ksummary \\
        -o ${prefix}_summaryline.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        Phoenix_summary_line.py: \$(${ica}Phoenix_summary_line.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
