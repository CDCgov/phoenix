def VERSION = '1.0' // Version information not provided by tool on CLI

process getMLST_scheme {
    tag "${meta.id}"
    label 'process_low'

    /*container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
        'quay.io/biocontainers/srst2:0.2.0--py27_2'}"
    */

    input:
    tuple val(meta)
    path(taxonomy)

    output:
    tuple val(meta), path("getmlst.out")               , emit: getMLST_out
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    test_title=\$(tail -n2 ${taxonomy} | head -n1 | cut -d\$'\t' -f1)
    if [[ \$test_title = "G:" ]]; then
        species=\$(tail -n1 ${taxonomy} | cut -d\$'\t' -f2)
        genus=\$(tail -n2 ${taxonomy} | head -n1 | cut -d\$'\t' -f2)
    elif [[ \$test_title = "s:" ]]; then
        species=\$(tail -n2 ${taxonomy} | head -n1 | cut -d\$'\t' -f2)
        genus=\$(tail -n3 ${taxonomy} | head -n1 | cut -d\$'\t' -f2)
    fi
    echo "\${genus} ___ \${species}"
    python -V

    /*  create lookup to insert here */
    db_entry="\${genus} \${species}"

    getMLST2.py --species '\$db_entry'

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getMLST: $VERSION
    END_VERSIONS
    """
}