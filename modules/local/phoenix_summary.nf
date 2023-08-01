process GATHER_SUMMARY_LINES {
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.0.0'

    input:
    path(summary_line_files)
    path(outdir_path)
    val(busco_val)

    output:
    path('Phoenix_Summary.tsv'), emit: summary_report
    path("versions.yml")             , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def busco_parameter = busco_val ? "--busco" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    Create_phoenix_summary_tsv.py \\
        --out Phoenix_Summary.tsv \\
        $busco_parameter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
