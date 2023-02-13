def VERSION = '1.1' // Version information not provided by tool on CLI

process CHECK_MLST {
    tag "${meta.id}"
    label 'process_low'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
    //    'quay.io/biocontainers/srst2:0.2.0--py27_2'}"

    container "quay.io/jvhagey/phoenix:base_v1.0.1"
    //container "quay.io/biocontainers/python:2.7--1"

    input:
    tuple val(meta), path(mlst_file), path(srst2_file), path(taxonomy_file), path(status)

    output:
    tuple val(meta), path("*_combined.tsv")                                                   , emit: checked_MLSTs
    tuple val(meta), path("*_status.txt")                                                     , emit: status
    path "versions.yml"                                                                       , emit: versions

    when:
    (task.ext.when == null || task.ext.when) && "${status[0]}" == "True"

    script:
    """
    wget $terra --secure-protocol=TLSv1_3 "https://pubmlst.org/data/dbases.xml"

    if [[ "${status[0]}" == "True" ]]; then
        check_and_fix_MLST2_new2.py --input $mlst_file --srst2 $srst2_file --taxonomy $taxonomy_file --docfile dbases.xml
    elif [[ "${status[0]}" == "False" ]]; then
        check_and_fix_MLST2_new2.py --input $mlst_file --taxonomy $taxonomy_file --docfile dbases.xml
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        check_mlst: $VERSION
    END_VERSIONS
    """
}
