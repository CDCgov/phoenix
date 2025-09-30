process ABRITAMR_REPORT {
    tag "$meta.id"
    label 'process_medium'
    container "quay.io/biocontainers/abritamr:1.0.17--pyh5707d69_1"

    input:
    tuple val(meta), path(summary_matches), path(summary_partials), val(organism_param)

    output:
    tuple val(meta), path("ID_${meta.id}_abritamr_summary.xlsx"), emit: summary
    tuple val(meta), path("abritamr.log")                       , emit: log
    tuple val(meta), path("mock_qc.txt")                        , emit: qc
    path "versions.yml"                          , emit: versions

    script:
    if ( "${organism_param[0]}" != "No Match Found") { organism_nf = "${organism_param[0]}" } else { organism_nf = "" }
    """
    # Assign the Nextflow variable to a Bash variable
    organism="$organism_nf"
    echo "ISOLATE,SPECIES_EXP,SPECIES_OBS,TEST_QC" > mock_qc.txt
    echo "ID_${meta.id},\${organism%%_*},\${organism%%_*},PASS" >> mock_qc.txt

    abritamr report \\
        --qc mock_qc.txt \\
        --runid ID_${meta.id}_abritamr_summary \\
        --matches ${summary_matches} \\
        --partials ${summary_partials}

    # fixing naming
    mv ID_${meta.id}_abritamr_summary_.xlsx ID_${meta.id}_abritamr_summary.xlsx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}