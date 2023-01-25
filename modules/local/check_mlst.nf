def VERSION = '1.0' // Version information not provided by tool on CLI

process CHECK_MLST {
    tag "${meta.id}"
    label 'process_low'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
    //    'quay.io/biocontainers/srst2:0.2.0--py27_2'}"
    //container "quay.io/biocontainers/python:2.7--1"

    container "quay.io/jvhagey/phoenix:base_v1.0.0" // This script uses python3. Python 3 is necessary to work in Terra

    input:
    tuple val(meta), path(mlst_file), path(taxonomy_file)

    output:
    tuple val(meta), path("*_combined.tsv")                                                   , emit: checked_MLSTs
    path "versions.yml"                                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_and_fix_MLST2.py --input $mlst_file --taxonomy $taxonomy_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        check_mlst: $VERSION
    END_VERSIONS
    """
}
