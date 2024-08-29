process CDIFF_RIBOTYPER {
    tag "$meta.id"
    label 'process_single'
    // v1.0.1 - MUST manually change below (line 27)!!!
    container 'quay.io/jvhagey/newtype@sha256:ca1edc26d6f64f62b3d4e09f167f15965d86520f630c23f0ceb6fd97fb3cab7e'

    input:
    tuple val(meta), path(csv_core)
    tuple val(meta), path(csv_accessory)

    output:
    tuple val(meta), path("*_ribotype.tsv"), emit: ribotype_file
    path("versions.yml")                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "v1.0.1"
    def container = task.container.toString() - "quay.io/jvhagey/newtype@"
    """
    # Call the real internal scripts to infer the ribotpes
    newtype.py -i ./ -s PN2.0 -o ${prefix}_ribotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        newtype.py: \$(newtype.py --version)
        container_tag: ${container_version}
        newtype_container: ${container}
    END_VERSIONS
    """
}
