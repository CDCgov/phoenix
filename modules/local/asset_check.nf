process ASSET_CHECK {
    label 'process_low'


    input:
    path(directory)


    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gunzip ${directory}/*.gz
    """
}