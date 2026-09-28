process BANDAGE {
    tag "${meta.id}"
    label 'process_medium'
    container 'staphb/bandage@sha256:3765e24bbbd7bdd1e34ef28825ca7303370d431901bba077f15cc07c7946a162'

    input:
    tuple val(meta), path(assembly_graph)

    output:
    tuple val(meta), path("*_bandage_graph.png"), emit: bandage_summary
    path ("versions.yml"),                        emit: versions

    script:
    def container = task.container.toString() - "staphb/bandage@"
    """
    Bandage image ${assembly_graph} ${meta.id}_bandage_graph.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bandage: \$(Bandage --version | sed -e "s/Version://g" )
        bandage_container: ${container}
    END_VERSIONS
    """

}