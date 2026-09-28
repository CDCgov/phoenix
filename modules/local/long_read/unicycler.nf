process UNICYCLER {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/unicycler@sha256:f611ddb4361f1151847de9d7e9b61ad23a5d197a1ec74c6237dc4b6107107a53'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}/assembly.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}/assembly.gfa"),   emit: gfa
    path "versions.yml",                                emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "staphb/unicycler@"
    """
    unicycler -l $fastq -o ${meta.id} -t $task.cpus --mode conservative

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        #unicycler: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//; s/ .*\$//')
        unicycler: \$( unicycler --version | cut -f 2 -d ' ' )
        unicycler_container: ${container}
    END_VERSIONS
    """
}
