process CDIFF_CLADE {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 23)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(mlst_combined_file), val(taxa_id), path(mlst_database)

    output:
    tuple val(meta), path("*clade.tsv"), emit: clade
    path "versions.yml"                , emit: versions

    when:
   "${taxa_id}" == "Clostridioides difficile"

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    def container_version = "base_v2.1.0"
    //def outdir_path = meta.project_id
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}get_Cdiff_clade.sh \\
        -m ${mlst_combined_file} \\
        -r "${mlst_database}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_Cdiff_clade.sh: \$(${ica}get_Cdiff_clade.sh -V)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
