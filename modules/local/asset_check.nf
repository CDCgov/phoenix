process ASSET_CHECK {
    label 'process_low'
    // base_v2.3.0 - MUST manually change below (line 22)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    path(zipped_sketch)
    path(mlst_db_path)
    path(kraken_db)
    path(clia_db_zipped)

    output:
    path('*.msh'),                          emit: mash_sketch
    path("versions.yml"),                   emit: versions
    path('db'),                             emit: mlst_db
    path('*_folder'),                       emit: kraken_db
    path('amrfinderdb_v*'), optional: true, emit: clia_db

    when:
    task.ext.when == null || task.ext.when

    script:
    def kraken_db_path = kraken_db ? "${kraken_db}" : "false" //checking if its null or an empty list
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def unzipped_sketch = "${zipped_sketch}".minus(".bz2")
    def unzip_clia_db = params.mode_upper == "CLIA" ? "tar --use-compress-program='pigz -vdf' -xf ${clia_db_zipped}" : "" 
    // Allow for multitude of zipped sources and remove the last extension, nevermind not needed for xz
    // def kraken_db_path = (!kraken_db || kraken_db.size() == 0) ? "false" : "${kraken_db}" //checking if its null or an empty list
    """
    ${unzip_clia_db}

    if [[ ${zipped_sketch} = *.gz ]]
    then
        echo "Unzipping gz file ${zipped_sketch}"
        pigz -vdf ${zipped_sketch}
        #for bz2 files
        #pigz -dc -L ${zipped_sketch} > ${unzipped_sketch}
    elif [[ ${zipped_sketch} = *.xz ]]
    then 
        echo "Unzipping xz file ${zipped_sketch}"
        xz -kfd --no-warn ${zipped_sketch}
    else
        :
    fi
    if [[ ${mlst_db_path} = *.tar.gz ]]
    then
        echo "Decompressing tar file ${mlst_db_path}"
        tar --use-compress-program="pigz -vdf" -xf ${mlst_db_path}
    else
        :
    fi
    if [[ ${kraken_db_path} != false ]]
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