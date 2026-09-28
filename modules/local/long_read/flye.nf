process FLYE {
    tag "${meta.id}"
    label 'process_high'
    container 'staphb/flye:2.9.3'

    input:
    tuple val(meta), path(fastq), path(fairy_outcome)

    output:
    tuple val(meta), path("*.fasta"),                 optional:true, emit: fasta
    tuple val(meta), path("*.txt"),                   optional:true, emit: assembly_info
    tuple val(meta), path("*.gfa"),                   optional:true, emit: gfa
    tuple val(meta), path("*.tsv"),                   optional:true, emit: flye_stat
    tuple val(meta), path("*.out"),                   optional:true, emit: flye_summary
    tuple val(meta), path("*_flye_outcome.csv"),                     emit: flye_outcome
    tuple val(meta), path('*_scaffolds_summary.txt'),                emit: outcome
    path ("versions.yml"),                                           emit: versions

    script:
    def container = task.container.toString() - "staphb/flye@"
    """
    # Set default that flye was successful
    # Lets downstream process know that flye completed ok - see assembly_failure.nf subworkflow
    flye_complete=run_failure
    echo \$flye_complete | tr -d "\\n" > ${meta.id}_flye_outcome.csv

    {

        # Try to assemble with flye. 
        flye --nano-hq $fastq -o . --iterations 1 --threads $task.cpus --scaffold || true

        if [[ -s assembly.fasta ]]; then

            # --- Flye actually produced an assembly: check whether scaffolds or contigs were made ---

            # Extract scaffold connection count
            SCAFFOLD_LINE=\$(grep -oE "Added [1-9]+ scaffold connections" flye.log | tail -n1 || true)

            # check if scaffolds or contigs were created and rename output files accordingly to stay inline with phx code
            # mirrors afterspades.sh: appends ,scaffolds_created/,no_scaffolds and ,contigs_created/,no_contigs
            if [[ -z "\$SCAFFOLD_LINE" ]]; then
                echo "Warning: could not find a 'scaffold connections' line in flye.log. Defaulting to treating the assembly as contigs-only" > flye.out
                mv assembly.fasta ${meta.id}.contigs.fasta
                # Compress to publish to outdir
                gzip ${meta.id}.contigs.fasta -c > ${meta.id}.contigs.fasta.gz

                echo ,no_scaffolds | tr -d "\\n" >> ${meta.id}_flye_outcome.csv
                echo ,contigs_created | tr -d "\\n" >> ${meta.id}_flye_outcome.csv
            else
                echo "Found log entry: '\$SCAFFOLD_LINE'" > flye.out
                mv assembly.fasta ${meta.id}.scaffolds.fasta
                gzip ${meta.id}.scaffolds.fasta -c > ${meta.id}.scaffolds.fasta.gz

                echo ,scaffolds_created | tr -d "\\n" >> ${meta.id}_flye_outcome.csv
                echo ,no_contigs | tr -d "\\n" >> ${meta.id}_flye_outcome.csv
            fi

            # rename output files to include the meta.id
            mv assembly_info.txt ${meta.id}.assembly_info.txt

            # Overwrite default that flye was successful to let downstream process know that flye completed ok - see assembly_failure.nf subworkflow
            flye_complete=run_completed
            echo \$flye_complete | tr -d "\\n" > ${meta.id}_flye_outcome.csv

        else

            # --- No assembly.fasta produced: Flye failed, regardless of its exit code ---

            echo "Flye pipeline ABORTED or FAILED -- no assembly.fasta was produced. See flye.log." > flye.out

            if [[ -f flye.log ]] && grep -q "No disjointigs were assembled" flye.log; then
                echo "Reason: No disjointigs were assembled -- likely insufficient reads/overlaps (too little data, wrong read type, or genome size mismatch)." >> flye.out
            elif [[ -f flye.log ]]; then
                echo "Error: pipeline aborted for an UNKNOWN reason. Check flye.log for the specific ERROR line(s):" >> flye.out
                grep "ERROR" flye.log >> flye.out || true
            else
                echo "Error: flye.log was not found -- flye may have crashed before writing any log output." >> flye.out
            fi

            flye_complete=run_failure
            echo \$flye_complete | tr -d "\\n" > ${meta.id}_flye_outcome.csv
            echo ,no_scaffolds | tr -d "\\n" >> ${meta.id}_flye_outcome.csv
            echo ,no_contigs | tr -d "\\n" >> ${meta.id}_flye_outcome.csv

        fi

    }

    if [[ -f ${meta.id}.assembly_info.txt ]]; then
        num_seqs=\$(tail -n +2 ${meta.id}.assembly_info.txt | grep -c '[^[:space:]]')
        if [[ "\$num_seqs" -gt 0 ]]; then
            printf "PASSED: More than 0 scaffolds in ${meta.id}." >> ${meta.id}_summary_old_4.txt
            cp ${meta.id}_summary_old_4.txt ${meta.id}_scaffolds_summary.txt
        else
            printf "FAILED: No scaffolds in ${meta.id}!" >> ${meta.id}_summary_old_4.txt
            cp ${meta.id}_summary_old_4.txt ${meta.id}_scaffolds_summary.txt
        fi
    else
        printf "FAILED: No scaffolds in ${meta.id}!" >> ${meta.id}_summary_old_4.txt
        cp ${meta.id}_summary_old_4.txt ${meta.id}_scaffolds_summary.txt
    fi


    # rename output files to include the meta.id, if it exists
    [[ -f flye.log ]] && mv flye.log ${meta.id}_flye.log

    # rename flye.out to include the meta.id, if it exists
    [[ -f flye.out ]] && mv flye.out ${meta.id}_flye.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version | sed -e "s/FLYE v//g" )
        flye_container: ${container}
    END_VERSIONS
    """
}