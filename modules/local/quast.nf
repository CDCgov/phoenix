process QUAST {
    tag "$meta.id"
    label 'process_low'
    container 'staphb/quast:5.0.2'

    input:
    tuple val(meta), path(consensus), val(fairy_outcome)

    output:
    tuple val(meta), path("quast")        , emit: results
    tuple val(meta), path('*.tsv')        , emit: report_tsv
    path "versions.yml"                   , emit: versions

    when:
    //if the files are not corrupt and there are equal number of reads in each file then run bbduk
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

    script:
    def args     = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    """
    quast.py \\
        --output-dir quast \\
        --threads $task.cpus \\
        $args \\
        $consensus

    mv quast/report.tsv ./${prefix}_summary.tsv

    # clean up name in file - allows multiQC to keep the sample together
    # Extract the current assembly name
    current_assembly=\$(awk '{print \$2}' "${prefix}_summary.tsv")

    # Check if the prefix matches the current assembly name
    if [[ "\$current_assembly" != *"${prefix}"* ]]; then
        # If not, update the file content
        sed -i "s/\\(Assembly\\s\\+\\).*.filtered.scaffolds/\\1${prefix}.filtered.scaffolds/" "${prefix}_summary.tsv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        quast: \$(quast.py --version 2>&1 | grep "QUAST" | sed 's/^.*QUAST v//; s/ .*\$//')
    END_VERSIONS
    """
}
