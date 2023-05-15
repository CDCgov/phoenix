def VERSION = '1.1' // Version information not provided by tool on CLI

process CHECK_MLST {
    tag "$meta.id"
    label 'process_single'
    container "quay.io/jvhagey/phoenix:base_v1.1.0"

    input:
    tuple val(meta), path(mlst_file), path(taxonomy_file)
    path(local_dbases)

    output:
    tuple val(meta), path("*_combined.tsv"), emit: checked_MLSTs
    tuple val(meta), path("*_status.txt"),   emit: status
    path("versions.yml")                 ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    fix_MLST2.py --input $mlst_file --taxonomy $taxonomy_file --mlst_database ${local_dbases}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        check_mlst: $VERSION
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
