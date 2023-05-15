process ASSET_CHECK {
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(zipped_sketch)
    path(mlst_db_path)

    output:
    path('*.msh'),        emit: mash_sketch
    path("versions.yml"), emit: versions
    path('db'),           emit: mlst_db

    when:
    task.ext.when == null || task.ext.when

    script:
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    if [[ ${zipped_sketch} = *.gz ]]
    then
        gunzip --force ${zipped_sketch}
    else
        :
    fi

    if [[ ${mlst_db_path} = *.gz ]]
    then
        tar -zvxf ${mlst_db_path}
    else
        :
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
