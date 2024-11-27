process BANDAGE {
    tag "$meta"
    label 'process_medium'
    errorStrategy 'ignore'
    container 'quay.io/biocontainers/bandage:0.8.1--hc9558a2_2'
    //sha256:79ae0ef6de06b68476667458da712bcddd2307bfd69d4cd3061944ab4c448ece

    input:
    tuple val(meta), path(gfa)

    output:
    tuple val(meta), path("*.png") , emit: bandage_summary
    path "versions.yml"            , emit: versions

    script:
    """
    Bandage image $gfa ${meta}_mqc.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    bandage: \$( bandage --version | sed -e "s/BANDAGE v//g" )
    END_VERSIONS
    """

}