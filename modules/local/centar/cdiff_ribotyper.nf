process CDIFF_RIBOTYPER {
    tag "$meta.id"
    label 'process_single'
    // v1.0.1 - MUST manually change below (line 27)!!!
    container 'quay.io/jvhagey/newtype@sha256:0b1b49acc633f88a071fb9157de405b88fb5e5aa4016e3afed46e43c4db03b30'

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
