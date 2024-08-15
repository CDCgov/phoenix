process COLLECT_PROJECT_FILES {
    tag "${meta.id}"
    stageInMode 'copy' // default is symlink. if its not set to copy changes in this script then changes original files.
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), file(software_versions), file(griphin_summary_tsv), file(griphin_summary_excel), file(phx_summary)
    //tuple val(meta), path(software_versions), path(griphin_summary_tsv), path(griphin_summary_excel), path(phx_summary)

    output:
    tuple val(meta), path("*_GRiPHin_Summary.tsv"),   emit: griphin_tsv
    tuple val(meta), path("*_GRiPHin_Summary.xlsx"),  emit: griphin_excel
    tuple val(meta), path("Phoenix_Summary.tsv"),     emit: phoenix_tsv
    tuple val(meta), path("software_versions.yml"),   emit: software_versions_file
    path("versions.yml"),                             emit: versions

    script: 
    // define variables
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    // Remove trailing slash if present
    //def clean_dir = dir.endsWith('/') ? "${dir}"[0..-2] : "${dir}"

    // Define project_path and ab_path
    //def project_path = clean_dir
    //def ab_path = "${clean_dir}/${meta.id}"
    """
    #just moving files so the output is accessible 
    mv input.*/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}