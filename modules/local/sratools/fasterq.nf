process SRATOOLS_FASTERQDUMP {
    tag "${meta.id}"
    label 'process_low'
    // 3.0.3--h87f3376_0 "quay.io/biocontainers/sra-tools@sha256:c9f92683e10091c3ef93066b1fcbdeeba89af49242ab778a9c8cc006f6be82a3"
    // 3.0.9--h9f5acd7_0
    container "quay.io/biocontainers/sra-tools@sha256:bd3dafdfb9ad5f301b72c5fdbbcbf411b19c59e117eabe4209dc15546c851c37"

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