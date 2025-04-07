process CDIFF_RIBOTYPER {
    tag "$meta.id"
    label 'process_medium'
    // v1.0.1 - MUST manually change below (line 27)!!!
    //container 'quay.io/jvhagey/newtype@sha256:c6e1aa3022330e0cf645e46523ea7de3706a04f03e5a3b24ffdb24f9d7b2d63c'
    container 'quay.io/enteevlachos/rosetta@sha256:6d4183671c71553010ac7ef4038b77402a05596b58f4c504a5b29da54cc35b59'

    input:
    tuple val(meta), path(csv_core)
    tuple val(meta), path(csv_accessory)

    output:
    tuple val(meta), path("*_ribotype.tsv"),               emit: ribotype_file
    tuple val(meta), path("*_ribotype_DetailedRport.tsv"), emit: detailed_ribotype_file
    path("versions.yml"),                                  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "v1.0.1"
    def container = task.container.toString() - "quay.io/jvhagey/newtype@"
    """
    mkdir ${prefix}
    mv *.gz ${prefix}

    # Call the real internal scripts to infer the ribotypes
    python3 /data/Rosetta_direct.py -i ${prefix} -s PN2.0 -o ${prefix}_ribotype.tsv

    newtype_version=\$(python3 /data/Rosetta_direct.py --version)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        newtype.py: \${newtype_version}
        container_tag: ${container_version}
        newtype_container: ${container}
    END_VERSIONS
    """
}
