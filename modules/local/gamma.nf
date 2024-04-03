process GAMMA {
    tag "$meta.id"
    label 'process_high'
    // v2.2 -- have to manually edit below (line 25)!!!
    container 'staphb/gamma@sha256:60d8ac58e016349a856fb7b443dd422ba69bae3f40e0dad83460d25ecf71101e'

    input:
    tuple val(meta), path(fasta), val(fairy_outcome)
    path(db)

    output:
    tuple val(meta), path("*.gamma")                , emit: gamma
    tuple val(meta), path("*.psl")                  , emit: psl
    tuple val(meta), path("*.gff")  , optional:true , emit: gff
    tuple val(meta), path("*.fasta"), optional:true , emit: fasta
    path "versions.yml"                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gamma_version = 2.2
    def container_version = task.container.toString() - "staphb/gamma@"
    """
    db_name=\$(echo $db | sed 's:.*/::' | sed 's/.fasta//')
    if [[ ${fasta} == *.gz ]]
    then
        FNAME=\$(basename ${fasta} .gz)
        gunzip -f ${fasta}
        GAMMA.py \\
        $args \\
        \$FNAME \\
        $db \\
        ${prefix}_\$db_name
    else
        GAMMA.py \\
        $args \\
        $fasta \\
        $db \\
        ${prefix}_\$db_name
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gamma: ${gamma_version}
        gamma_container: ${container_version}
        Database: $db
    END_VERSIONS
    """
}
