process COMBINE_SRST2_MLST {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta), path(line_files)   // collected list of *_srst2_line.txt (may be empty list if all no-match)

    output:
    tuple val(meta), path("*_srst2_temp.mlst"), emit: mlst_results_temp
    tuple val(meta), path("*_srst2_status.txt"), emit: empty_checker

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def header = "Sample\tdatabase\tST\tmismatches\tuncertainty\tdepth\tmaxMAF\tlocus_1\tlocus_2\tlocus_3\tlocus_4\tlocus_5\tlocus_6\tlocus_7\tlocus_8\tlocus_9\tlocus_10"
    """
    echo -e "${header}" > ${prefix}_srst2_temp.mlst

    # line_files may be a single path or a list; Nextflow stages them all into cwd either way.
    # Concatenate in filename order for determinism (sorted by scheme name), since arrival
    # order is not guaranteed under parallel execution and row order was never semantically
    # meaningful — but a deterministic order still aids reproducibility/debugging.
    for f in \$(ls *_srst2_line.txt 2>/dev/null | sort); do
        cat "\${f}" >> ${prefix}_srst2_temp.mlst
    done

    if [[ \$(wc -l < ${prefix}_srst2_temp.mlst) -gt 1 ]]; then
        echo "True" > ${prefix}_srst2_status.txt
    else
        echo "False" > ${prefix}_srst2_status.txt
    fi
    """
}