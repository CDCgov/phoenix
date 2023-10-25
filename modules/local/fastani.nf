process FASTANI {
    tag "$meta.id"
    label 'process_medium'
    container 'staphb/fastani:1.33'
    //stageInMode "copy"
    /*container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastani:1.33--h0fdf51a_0' :
        'quay.io/biocontainers/fastani:1.33--h0fdf51a_0' }"*/

    input:
    tuple val(meta), path(query), path(reference), path(reference_dir)

    output:
    tuple val(meta), path("*.ani.txt"), emit: ani
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    line=\$(head -n1 ${reference})
    if [[ "\${line}" = "No MASH hit found" ]]; then
        echo "Mash/FastANI Error: No MASH hit found" > ${prefix}.ani.txt
    else
        db_version=\$(echo ${reference} | sed 's/_best_MASH_hits.txt//' | sed 's/${prefix}_//' )
        # Setup to catch any issues while grabbing date from DB name
        if [[ "\${db_version}" = "" ]]; then
            db_version="REFSEQ_unknown"
        fi
        fastANI \\
            -q $query \\
            --rl $reference \\
            -o ${prefix}_\${db_version}.ani.txt
    fi

    # NOTE: No ANI output is reported (but file is created) for a genome pair if ANI value is much below 80%. 
    # Such case should be computed at amino acid level. However, we aren't going to do that so we will confirm the file isn't empty. 
    # if the file doesn't exist or is empty say no Mash hits found
    if [[ ! -s ${prefix}_\${db_version}.ani.txt ]]; then
        echo "Mash/FastANI Error: No hits above an ANI value >80%" > ${prefix}_\${db_version}.ani.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani: \$(fastANI --version 2>&1 | sed 's/version//;')
    END_VERSIONS
    """
}
