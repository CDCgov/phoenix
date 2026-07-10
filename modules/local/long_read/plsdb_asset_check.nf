process PLSDB_ASSET_CHECK {
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 22)!!!
    //container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'
    container 'ncbi/blast@sha256:81f118d2e4f7e11494d27fdbb99c9430423105afff50c4ae158db41d58a3fc57'  //ncbi/blast:2.17.0

    input:
    // path(zipped_sketch)
    // path(mlst_db_path)
    // path(kraken_db)
    // path(clia_db_zipped)
    //path(zipped_plsdb)
    path(plsdb_dir)

    output:
    // path('*.msh'),                          emit: mash_sketch
    path("versions.yml"),                   emit: versions
    // path('db'),                             emit: mlst_db
    // path('*_folder'),                       emit: kraken_db
    // path('amrfinderdb_v*'), optional: true, emit: clia_db
    path("${plsdb_dir}/plsdb_2024_05_31_v2.fasta"),                           emit: plsdb

    // when:
    // task.ext.when == null || task.ext.when 

    script:
    // def kraken_db_path = kraken_db ? "${kraken_db}" : "false" //checking if its null or an empty list
    // def container_version = "base_v2.2.0"
    // def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    // def unzipped_sketch = "${zipped_sketch}".minus(".bz2")
    // def unzip_clia_db = params.mode_upper == "CLIA" ? "tar --use-compress-program='pigz -vdf' -xf ${clia_db_zipped}" : "" 
    def blast_version = "v2.17.0"
    def blast_container = task.container.toString() - "ncbi/blast@"
    //def unzipped_plsdb = "${zipped_plsdb}".minus(".gz")

    // Allow for multitude of zipped sources and remove the last extension, nevermind not needed for xz
    // def kraken_db_path = (!kraken_db || kraken_db.size() == 0) ? "false" : "${kraken_db}" //checking if its null or an empty list
    """
    
    if [[ ! -f "${plsdb_dir}/plsdb_2024_05_31_v2.fasta" && ! -f "${plsdb_dir}/plsdb_2024_05_31_v2.nsq" ]]
    then
        echo "Download plsdb v.2024_05_31_v2 from https://ccb-microbe.cs.uni-saarland.de/plsdb2025, unzip and make blast database"
        wget "https://ndownloader.figshare.com/files/49859772" --no-check-certificate -O ${plsdb_dir}/plsdb_2024_05_31_v2.fasta.bz2
        bzip2 -d ${plsdb_dir}/plsdb_2024_05_31_v2.fasta.bz2
        makeblastdb -in ${plsdb_dir}/plsdb_2024_05_31_v2.fasta -parse_seqids -dbtype nucl -title "PLSDB 2024_05_31_v2" -out ${plsdb_dir}/plsdb_2024_05_31_v2
    else
        :
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast_version: ${blast_version}
        blast_container: ${blast_container}
    END_VERSIONS
    """
}
