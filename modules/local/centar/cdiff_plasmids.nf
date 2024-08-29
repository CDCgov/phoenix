process CDIFF_PLASMIDS {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 27)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(assembly), val(fairy_outcome)
    path(plasmid_database)

    output:
    tuple val(meta), path("*.tsv"), emit: plasmids_file
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "jvhagey/phoenix@"
    """
    # Call the real internal scripts to infer the ribotpes

    # Placeholder to create a blank file  
    touch "${prefix}_plasmids.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plasmid script: 0.0.0
        plasmid_container: ${container}
    END_VERSIONS
    """
}
