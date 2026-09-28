process SEQKIT_RAWSTATS {
    tag "${meta.id}"
    label 'process_medium'
    container 'staphb/seqkit@sha256:1841442f1bc6bfd0a711e5a1465b1e350ad7a57f037634bbddb4b003381da410'

    input:
    tuple val(meta), path(reads), path(fairy_outcome)

    output:
    tuple val(meta), path("*_LR_raw_read_counts.txt"), emit: rawstats
    tuple val(meta), path('*_summary_old_2.txt'),      emit: outcome_to_edit
    path ("versions.yml"),                             emit: versions

    script:
    def container = task.container.toString() - "staphb/seqkit@"
    """
    echo "PASSED: Long-read no read pairs." >> ${meta.id}_summary_old.txt

    seqkit stats --tabular --threads $task.cpus --all ${reads} > ${meta.id}_LR_raw_read_counts.txt

    # making a copy of the summary file to pass downstream to handle file names being the same
    cp ${meta.id}_summary_old.txt ${meta.id}_summary_old_2.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed -e "s/seqkit //g")
        seqkit_container: ${container}
    END_VERSIONS
    """
}
