process QUAST {
    tag "$meta.id"
    label 'process_medium'
    //5.3.0
    //container 'staphb/quast:5.0.2'
    container 'staphb/quast@sha256:850dd9821c0c92c8bb1b258658a64b39682e8f3d61dd1fa6f40e3ef906f1edb8'

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
    //set up for terra
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = "PATH=/opt/conda/envs/quast/bin:\$PATH"
        terra_exit = """PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/quast/bin:||')" """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
    #adding path for running quast on terra
    $terra

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

    #revert path back to main envs for running on terra
    $terra_exit
    """
}
