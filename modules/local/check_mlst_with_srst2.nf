process CHECK_MLST_WITH_SRST2 {
    tag "$meta.id"
    label 'process_single'
    container "quay.io/jvhagey/phoenix:base_v2.0.2"

    input:
    tuple val(meta), path(mlst_file), path(srst2_file), path(taxonomy_file), val(status), path(local_dbases)

    output:
    tuple val(meta), path("*_combined.tsv"), emit: checked_MLSTs
    tuple val(meta), path("*_status.txt"),   emit: status
    path("versions.yml"),                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    if [[ "${status[0]}" == "True" ]]; then
        fix_MLST2.py --input $mlst_file --srst2 $srst2_file --taxonomy $taxonomy_file --mlst_database $local_dbases
    elif [[ "${status[0]}" == "False" ]]; then
        fix_MLST2.py --input $mlst_file --taxonomy $taxonomy_file --mlst_database $local_dbases
    else 
        echo "Something went very wrong, please open an issue on Github for the PHoeNIx developers to address."
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}