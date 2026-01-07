process ABRITAMR {
    tag "$meta.id"
    label 'process_medium'
    container "quay.io/biocontainers/abritamr:1.0.17--pyh5707d69_1"

    input:
    tuple val(meta), path(fasta), val(organism_param), path(ar_bd)

    output:
    tuple val(meta), path("*.summary_matches.txt")  , emit: matches
    tuple val(meta), path("*.summary_partials.txt") , emit: partials
    tuple val(meta), path("*.summary_virulence.txt"), emit: virulence
    tuple val(meta), path("*.amrfinder.out")        , emit: amrfinder
    tuple val(meta), path("*.abritamr.txt")         , emit: summary, optional: true
    path "versions.yml"                             , emit: versions

    script:
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    fasta_name = fasta.getName().replace(".gz", "")
    // use --species
    if ( "${organism_param[0]}" != "No Match Found") { species = "--species ${organism_param[0]}" } else { species = "" }
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    #have to add ID in front for handling for if the ID is only numbers
    abritamr run \\
        --contigs $fasta_name \\
        --prefix ${meta.id} \\
        --jobs $task.cpus \\
        --amrfinder_db ${ar_bd} \\
        $species

    # Rename output files to prevent name collisions
    mv ${meta.id}/summary_matches.txt ./${meta.id}.summary_matches.txt
    mv ${meta.id}/summary_partials.txt ./${meta.id}.summary_partials.txt
    mv ${meta.id}/summary_virulence.txt ./${meta.id}.summary_virulence.txt
    mv ${meta.id}/amrfinder.out ./${meta.id}.amrfinder.out
    if [ -f ${meta.id}/abritamr.txt ]; then
        # This file is not always present
        mv ${meta.id}/abritamr.txt ./${meta.id}.abritamr.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' | sed 's/)//'))
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}