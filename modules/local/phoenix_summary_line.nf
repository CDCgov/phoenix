process CREATE_SUMMARY_LINE {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

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
    path('*_summaryline.tsv'), emit: line_summary
    path("versions.yml")     , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    // allowing for some optional parameters for -entry SCAFFOLDS/CDC_SCAFFOLDS nothing should be passed.
    def trimmed_qc_data = trimmed_qc_data_file ? "-t $trimmed_qc_data_file" : ""
    def trim_ksummary   = trimd_ksummary ? "-k $trimd_ksummary" : ""
    def fastani_file    = fastani ? "-f $fastani" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    Phoenix_summary_line.py \\
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
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
