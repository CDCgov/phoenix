process BANDAGE {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/bandage:0.8.1--hc9558a2_2'
    //sha256:79ae0ef6de06b68476667458da712bcddd2307bfd69d4cd3061944ab4c448ece

    input:
    tuple val(meta), path(assembly_graph)

    output:
    tuple val(meta), path("*_bandage_graph.png"), emit: bandage_summary
    path ("versions.yml"),                        emit: versions

    script:
    """
    Bandage image ${assembly_graph} ${meta.id}_bandage_graph.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bandage: \$(Bandage --version | sed -e "s/Version://g" )
    END_VERSIONS
    """

}