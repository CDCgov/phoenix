process CDIFF_PLASMIDS {
    tag "$meta.id"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 24)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

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
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
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
