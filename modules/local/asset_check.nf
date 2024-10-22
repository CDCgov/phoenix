process ASSET_CHECK {
    label 'process_low'
    // base_v2.1.0 - MUST manually change below (line 21)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    path(zipped_sketch)
    path(mlst_db_path)
    path(kraken_db)

    output:
    path('*.msh'),        emit: mash_sketch
    path("versions.yml"), emit: versions
    path('db'),           emit: mlst_db
    path('*_folder'),     emit: kraken_db

    when:
    task.ext.when == null || task.ext.when

    script:
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def unzipped_sketch = "${zipped_sketch}".minus(".bz2")
    def kraken_db_path = (!kraken_db || kraken_db.size() == 0) ? "false" : "${kraken_db}" //checking if its null or an empty list
    """
    if [[ ${zipped_sketch} = *.gz ]]
    then
        pigz -vdf ${zipped_sketch}
        #for bz2 files
        #pigz -dc -L ${zipped_sketch} > ${unzipped_sketch}
    else
        :
    fi

    if [[ ${mlst_db_path} = *.tar.gz ]]
    then
        tar --use-compress-program="pigz -vdf" -xf ${mlst_db_path}
    else
        :
    fi

    if [[ ${kraken_db_path} != "false" ]]
    then
        if [[ ${kraken_db_path} = *.tar.gz ]]
        then
            folder_name=\$(basename ${kraken_db_path} .tar.gz)
            tar --use-compress-program="pigz -vdf" -xf ${kraken_db_path}
            mkdir \${folder_name}_folder
            mv *.kmer_distrib \${folder_name}_folder
            mv *.k2d \${folder_name}_folder
            mv seqid2taxid.map \${folder_name}_folder
            mv inspect.txt \${folder_name}_folder
            mv ktaxonomy.tsv \${folder_name}_folder
        else
            folder_name=\$(basename ${kraken_db_path} .tar.gz)
            mv \${folder_name} \${folder_name}_folder
        fi
    else
        #just make an empty folder to keep things moving
        mkdir empty_folder
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}