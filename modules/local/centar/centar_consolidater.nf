process CENTAR_CONSOLIDATER {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 23)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(tox_file), path(clade_file)

    output:
    tuple val(meta), path("*clade.tsv"), emit: clade
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.terra==false) { terra = ""} 
    else if (params.terra==true) { terra = "-t terra" }
    else { error "Please set params.terra to either \"true\" or \"false\"" }
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Centar_Consolidater.sh \\
        -t ${tox_file} \\
        -c "${clade_file}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Centar_Consolidater.sh: \$(${ica}get_cdiff_clade.sh -V)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
