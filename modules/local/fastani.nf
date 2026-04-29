process FASTANI {
    tag "$meta.id"
    label 'process_medium'
    // v1.34
    container 'staphb/fastani@sha256:d2c1bf7f1792c7d371904e8744f1a974535f0814d29212b254411e78c1436f59'

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
    def container = task.container.toString() - "staphb/fastani@"
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

        # NOTE: No ANI output is reported (but file is created) for a genome pair if ANI value is much below 80%. 
        # Such case should be computed at amino acid level. However, we aren't going to do that so we will confirm the file isn't empty. 
        # if the file doesn't exist or is empty say no Mash hits found
        if [[ ! -s ${prefix}_\${db_version}.ani.txt ]]; then
            echo "Mash/FastANI Error: No hits above an ANI value >80%" > ${prefix}_\${db_version}.ani.txt
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastani:\$(fastANI --version 2>&1 | sed 's/version//;')
        fastani_container: ${container}
    END_VERSIONS
    """
}
