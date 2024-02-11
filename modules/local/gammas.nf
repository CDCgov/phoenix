process GAMMA_S {
    tag "$meta.id"
    label 'process_high'
    container 'staphb/gamma:2.2'
    /*container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gamma%3A2.1--hdfd78af_0':
        'quay.io/biocontainers/gamma:2.1--hdfd78af_0' }"*/

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path("*.gamma")                , emit: gamma
    tuple val(meta), path("*.psl")                  , emit: psl
    tuple val(meta), path("*.gff")  , optional:true , emit: gff
    tuple val(meta), path("*.fasta"), optional:true , emit: fasta
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = task.container.toString() - "staphb/gamma:"
    """
    db_name=\$(echo $db | sed 's:.*/::' | sed 's/.fasta//')
    if [[ ${fasta} == *.gz ]]
    then
        FNAME=\$(basename ${fasta} .gz)
        gunzip -f ${fasta}
        GAMMA-S.py \\
        $args \\
        \$FNAME \\
        $db \\
        ${prefix}_\$db_name
    else
        GAMMA-S.py \\
        $args \\
        $fasta \\
        $db \\
        ${prefix}_\$db_name
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gamma: ${container_version}
        Database: $db
    END_VERSIONS
    """
}
