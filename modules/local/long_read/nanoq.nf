process NANOQ {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/biocontainers/nanoq@sha256:86e6f65c7a0c8e626511c7608de484c0bcfa418f3bd39af12250fe61c02b0fdb'

    input:
    tuple val(meta), path(subfastq), path(fairy_outcome)
    val(length)
    val(qscore)

    output:
    tuple val(meta), path("*_trim.fastq.gz"),                        emit: fastq
    tuple val(meta), path("*_LR_trimmed_read_counts.txt"),           emit: trimmed_stats
    tuple val(meta), path('*_trimstats_summary.txt'), optional:true, emit: outcome
    tuple val(meta), path('*_summary_old_4.txt'),                    emit: outcome_to_edit
    path ("versions.yml"),                                           emit: versions

    script:
    def container = task.container.toString() - "quay.io/biocontainers/nanoq@"
    """
    nanoq --input $subfastq --min-len $length --min-qual $qscore --report ${meta.id}_LR_trimmed_read_counts.txt --stats --header --output-type g -o ${meta.id}_trim.fastq.gz

    #check if reads remaining after trimming, if not then report failure
    if [[ \$(awk 'NR==2{print \$1}' ${meta.id}_LR_trimmed_read_counts.txt) -eq 0 ]]; then
        echo "FAILED: There are 0 reads in ${meta.id} after trimming!" >> ${meta.id}_summary_old_3.txt
        cp ${meta.id}_summary_old_3.txt ${meta.id}_summary_old_4.txt
        cp ${meta.id}_summary_old_3.txt ${meta.id}_trimstats_summary.txt
    else
        echo "PASSED: There are reads in ${meta.id} after trimming!" >> ${meta.id}_summary_old_3.txt
        cp ${meta.id}_summary_old_3.txt ${meta.id}_summary_old_4.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoq: \$(nanoq --version | sed -e "s/nanoq //g")
        nanoq_container: ${container}
    END_VERSIONS
    """
}
