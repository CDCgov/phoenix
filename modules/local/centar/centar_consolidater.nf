process CENTAR_CONSOLIDATER {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 23)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    tuple val(meta), path(tox_file), 
    path(clade_file), 
    path(toxinotype_file), 
    path(other_AR_AA_file),
    path(other_AR_NT_file), 
    path(plasmids_file), 
    path(rt_file)
    path(st_rt_xwalk)

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
    def container_version = "base_v2..0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def ribotype_file = rt_file ? "-r ${rt_file}" : ""
    """
    ${ica}Centar_Consolidater.sh \\
        -t ${tox_file} \\
        -c ${clade_file} \\
        -y ${toxinotype_file} \\
        -a ${other_AR_AA_file} \\
        -n ${other_AR_NT_file} \\
        ${ribotype_file} \\
        -p ${plasmids_file} \\
        -o ${prefix}_centar_output.tsv \\
        -s ${prefix} \\
        -x ${st_rt_xwalk}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Centar_Consolidater.sh: \$(${ica}Centar_Consolidater.sh -V)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
