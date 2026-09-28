process LRGE {
    tag "${meta.id}"
    label 'process_medium'
    container 'staphb/lrge@sha256:7a74eca61650b4f50d01072bc4e61208bf4c87f9aa438879693072aec8b70bd4'

    input:
    tuple val(meta), path(fastq_lr), path(fairy_outcome)

    output:
    tuple val(meta), path("*_genome_size.txt"),                       emit: estimation
    tuple val(meta), path('*_estimation_summary.txt'), optional:true, emit: outcome
    tuple val(meta), path('*_summary_old_3.txt'),                     emit: outcome_to_edit
    path ("versions.yml"),                                            emit: versions

    script:
    def container = task.container.toString() - "staphb/lrge@"
    """
    # Run lrge without letting a non-zero exit kill the script outright. We need to inspect stderr first to decide error handling.
    set +e
    lrge -t $task.cpus ${fastq_lr} > ${meta.id}_genome_size.txt 2> lrge_stderr.txt
    lrge_exit_code=\$?
    set -e

    if [ \$lrge_exit_code -ne 0 ]; then
        if grep -q "No finite estimates were generated" lrge_stderr.txt; then
            # Expected/known failure mode -- lrge couldn't produce an estimate for this sample (e.g. insufficient read overlap). 
            # Leave the genome size file blank, log a line to the fairy outcome file, and let the pipeline continue.
            echo "FAILED: LRGE could not generate a genome size estimate for ${meta.id} (no finite estimates were generated)." >> ${meta.id}_summary_old_2.txt
        else
            # Any other failure reason is unexpected -- fail the pipeline as normal so it isn't silently erroring.
            echo "FAILED: LRGE failed for ${meta.id} for a reason other than 'No finite estimates were generated'" >> ${meta.id}_summary_old_2.txt
            exit \$lrge_exit_code
        fi
        # making a copy of the summary file to pass downstream to handle file names being the same
        cp ${meta.id}_summary_old_2.txt ${meta.id}_summary_old_3.txt
        cp ${meta.id}_summary_old_3.txt ${meta.id}_estimation_summary.txt
    else
        # LRGE succeeded -- log a line to the fairy outcome file.
        echo "PASSED: LRGE genome size estimate generated for ${meta.id}." >> ${meta.id}_summary_old_2.txt
        # making a copy of the summary file to pass downstream to handle file names being the same
        cp ${meta.id}_summary_old_2.txt ${meta.id}_summary_old_3.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LRGE: \$(LRGE --version | sed -e "s/LRGE//g")
        lrge_container: ${container}
    END_VERSIONS
    """
}
