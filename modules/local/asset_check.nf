process ASSET_CHECK {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(zipped_sketch)

    output:
    path('*.msh'), emit: mash_sketch

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [[ ${zipped_sketch} = *.gz ]]
    then
        gunzip --force ${zipped_sketch}
    else
        :
    fi
    """
}
