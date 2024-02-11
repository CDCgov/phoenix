process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    // v0.12.1
    container 'staphb/fastqc@sha256:f5d8f72753269e0cee071fe198c89a59a1f8071445739b3398f7818f7cb039ae'

    input:
    tuple val(meta), path(reads), val(fairy_outcome)

    when:
    //if there are scaffolds left after filtering
    "${fairy_outcome[3]}" == "PASSED: There are reads in ${meta.id} R1/R2 after trimming."

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path("versions.yml")           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def container = task.container.toString() - "staphb/fastqc@"
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        fastqc $args --threads $task.cpus ${prefix}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        fastqc $args --threads $task.cpus ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
            fastqc_container: ${container}
        END_VERSIONS
        """
    }
}
