process UNICYCLER {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/unicycler:0.5.0'
    //sha256:f1e556959e2b6df92d66726ed9743bb17fccd1a6a2bf4961dcc399512b03512d

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}/assembly.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}/assembly.gfa"),   emit: gfa
    path "versions.yml",                                emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    unicycler -l $fastq -o ${meta.id} -t $task.cpus --mode conservative

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        #unicycler: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//; s/ .*\$//')
        unicycler: \$( unicycler --version | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}
