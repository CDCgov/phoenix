process UPDATE_SUMMARY_FILE {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.3.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(phoenix_line_summary), \
    path(mlst_file), \
    path(ar_gamma_file), \
    path(amr_report),\
    path(synopsis)

    output:
    path('Updated_*_summaryline.tsv'), emit: updated_line_summary
    path('Updated_*_combined.tsv'),    emit: updated_mlst
    path("versions.yml")             , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    // allowing for some optional parameters for -entry SCAFFOLDS/CDC_SCAFFOLDS nothing should be passed.
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Update_Phoenix_summary_line.py \\
        -g ${ar_gamma_file} \\
        -a ${amr_report} \\
        -m ${mlst_file} \\
        -s ${phoenix_line_summary} \\
        -o Updated_${prefix}_summaryline.tsv

    # Extract the first element in a tab-separated line using cut
    Database=\$(sed -n '2p' "${mlst_file}" | cut -f4 | awk '{print toupper(\$0)}')
    ST=\$(sed -n '2p' "${mlst_file}" | cut -f5)
    Source=\$(sed -n '2p' "${mlst_file}" | cut -f2)
    new_mlst_line="MLST-\$database                 : SUCCESS  : ST\$ST via \$Source"

    # Use sed to replace the line starting with "MLST" with the new content
    sed -i '/^MLST/c\'"\$new_mlst_line" ${synopsis}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        Update_Phoenix_summary_line.py: \$(${ica}Update_Phoenix_summary_line.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
