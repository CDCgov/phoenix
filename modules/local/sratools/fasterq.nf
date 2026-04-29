process SRATOOLS_FASTERQDUMP {
    tag "${meta.id}"
    label 'process_low'
    // 3.2.0--h4304569_0
    container "quay.io/biocontainers/sra-tools@sha256:db636fa5785c482fe69d836b2c4e24c9a912b9557ed069cad3835d4234b9354e"

    input:
    tuple val(meta), path(sra_folder)

    output:
    tuple val(meta), path("*_*.fastq.gz"), emit: reads // we don't want the SRR.fastq just the forward and reverse
    path("versions.yml"),                  emit: versions

    script:
    //define variables
    def args = task.ext.args ?: ''
    def srr_number = sra_folder.toString() - "_Folder"
    def container = task.container.toString() - "quay.io/biocontainers/sra-tools@"
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
        sratools_container: ${container}
    END_VERSIONS
    """
}