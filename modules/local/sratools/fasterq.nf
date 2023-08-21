process SRATOOLS_FASTERQDUMP {
    tag "${meta.id}"
    label 'process_low'
    container "quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0"

    input:
    tuple val(meta), path(sra_folder)

    output:
    tuple val(meta), path("*_*.fastq.gz"), emit: reads // we don't want the SRR.fastq just the forward and reverse
    path("versions.yml"),                  emit: versions

    script:
    def args = task.ext.args ?: ''
    def srr_number = sra_folder.toString() - "_Folder"
    """
    # change folder name back for fasterq-dump to find
    mv ${sra_folder} ${srr_number}

    fasterq-dump \\
        $args \\
        --threads $task.cpus \\
        ${srr_number}

    gzip ${srr_number}_1.fastq
    gzip ${srr_number}_2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | sed 's/fasterq-dump : //' | awk 'NF' )
    END_VERSIONS
    """
}