process ABRITAMR {
    tag "$meta.id"
    label 'process_medium'
    container "quay.io/biocontainers/abritamr:1.0.17--pyh5707d69_1"

    input:
    tuple val(meta), path(fasta), val(fairy_outcome)

    output:
    tuple val(meta), path("*.summary_matches.txt")  , emit: matches
    tuple val(meta), path("*.summary_partials.txt") , emit: partials
    tuple val(meta), path("*.summary_virulence.txt"), emit: virulence
    tuple val(meta), path("*.amrfinder.out")        , emit: amrfinder
    tuple val(meta), path("*.abritamr.txt")         , emit: summary, optional: true
    path "*.{log,err}"                              , emit: logs, optional: true
    path ".command.*"                               , emit: nf_logs
    path "versions.yml"                             , emit: versions

    when:
    //if the files are not corrupt and there are equal number of reads in each file then run bbduk
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

    script:
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    abritamr run \\
        --contigs $fasta_name \\
        --prefix results \\
        --jobs $task.cpus

    # Rename output files to prevent name collisions
    mv results/summary_matches.txt ./${meta.id}.summary_matches.txt
    mv results/summary_partials.txt ./${meta.id}.summary_partials.txt
    mv results/summary_virulence.txt ./${meta.id}.summary_virulence.txt
    mv results/amrfinder.out ./${meta.id}.amrfinder.out
    if [ -f results/abritamr.txt ]; then
        # This file is not always present
        mv results/abritamr.txt ./${meta.id}.abritamr.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}