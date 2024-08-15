process CDIFF_TOXINOTYPER {
    tag "${meta.id}"
    label 'process_single'
    container 'staphb/gamma@sha256:60d8ac58e016349a856fb7b443dd422ba69bae3f40e0dad83460d25ecf71101e'

    input:
    tuple val(meta), path(assembly), val(fairy_outcome)
    path(tox_database)
    path(tox_definitions)

    output:
    tuple val(meta), path("*.tox"), emit: tox_file
    tuple val(meta), path("*.psl"), emit: tox_psl_file
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "staphb/gamma@"
    """
    ${ica}blat_toxinotypes.sh \\
        -i ${assembly} \\
        -d ${tox_database} \\
        -t ${tox_definitions} \\
        -o ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blat_toxinotypes.sh: \$(${ica}blat_toxinotypes.sh -V)
        gamma_container: ${container}
    END_VERSIONS
    """
}
