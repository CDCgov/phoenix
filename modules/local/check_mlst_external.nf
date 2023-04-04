def VERSION = '1.1' // Version information not provided by tool on CLI

process CHECK_MLST {
    tag "${meta.id}"
    label 'process_single'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
    //    'quay.io/biocontainers/srst2:0.2.0--py27_2'}"

    //container "quay.io/jvhagey/phoenix:base_v1.0.1"
    container "quay.io/biocontainers/python:2.7--1"

    input:
    tuple val(meta), path(mlst_file), path(srst2_file), path(taxonomy_file)

    output:
    tuple val(meta), path("*_combined.tsv")                                                   , emit: checked_MLSTs
    tuple val(meta), path("*_status.txt")                                                     , emit: status
    path "versions.yml"                                                                       , emit: versions

    when:
    (task.ext.when == null || task.ext.when)

    script:
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) {
        terra = ""
    } else if (params.terra==true) {
        terra = "--no-check-certificate"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    check_and_fix_MLST2_new2.py --input $mlst_file --srst2 $srst2_file --taxonomy $taxonomy_file
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        check_mlst: $VERSION
    END_VERSIONS
    """
}
