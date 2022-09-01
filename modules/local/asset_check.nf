process ASSET_CHECK {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.0.0'

    input:
    path(zipped_sketch)

    output:
    path('REFSEQ_20210820_Bacteria_complete.msh'), emit: mash_sketch

    when:
    task.ext.when == null || task.ext.when

    shell:
    """
    if [[ $zipped_sketch == *.gz ]]
    then
        gunzip --force $zipped_sketch
    else
        :
    fi
    """
}