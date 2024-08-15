process CENTAR_CONSOLIDATER {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 23)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(tox_file), 
    path(clade_file), 
    path(toxinotype_file), 
    path(other_AR_file), 
    path(rt_file), 
    path(plasmids_file)

    output:
    tuple val(meta), path("*.tsv"), emit: centar_summary_line
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Centar_Consolidater.sh \\
        -t ${tox_file} \\
        -c ${clade_file} \\
        -t ${toxinotype_file} \\
        -a ${other_AR_file} \\
        -r ${rt_file} \\
        -p ${plasmids_file}
        -o ${prefix}_centar_output.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(${ica}Centar_Consolidator.sh -V)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
