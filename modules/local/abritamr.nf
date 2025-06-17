process ABRITAMR {
    tag "$meta.id"
    label 'process_medium'
    container "quay.io/biocontainers/abritamr:1.0.17--pyh5707d69_1"

    input:
    tuple val(meta), path(fasta), val(fairy_outcome), val(organism_param)
    path(ar_bd)

    output:
    tuple val(meta), path("*.summary_matches.txt")  , emit: matches
    tuple val(meta), path("*.summary_partials.txt") , emit: partials
    tuple val(meta), path("*.summary_virulence.txt"), emit: virulence
    tuple val(meta), path("*.amrfinder.out")        , emit: amrfinder
    tuple val(meta), path("*.abritamr.txt")         , emit: summary, optional: true
    path "versions.yml"                             , emit: versions

    when:
    //if the files are not corrupt and there are equal number of reads in each file then run bbduk
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

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
        --prefix ID_${meta.id} \\
        --jobs $task.cpus \\
        --d ${ar_bd} \\
        $species

    # Rename output files to prevent name collisions
    mv ID_${meta.id}/summary_matches.txt ./ID_${meta.id}.summary_matches.txt
    mv ID_${meta.id}/summary_partials.txt ./ID_${meta.id}.summary_partials.txt
    mv ID_${meta.id}/summary_virulence.txt ./ID_${meta.id}.summary_virulence.txt
    mv ID_${meta.id}/amrfinder.out ./ID_${meta.id}.amrfinder.out
    if [ -f ${meta.id}/abritamr.txt ]; then
        # This file is not always present
        mv ID_${meta.id}/abritamr.txt ./ID_${meta.id}.abritamr.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}