process KRAKEN_DB_DWNLD {

    output:
    tuple val(meta), path('taxo.k2d')       , optional:false, emit: scaffolds
    tuple val(meta), path('hash.k2d')       , optional:false, emit: contigs
    tuple val(meta), path('opts.k2d')       , optional:false, emit: transcripts

    script:
    """
    cd ${baseDir}/assets/databases
    wget -c https://refdb.s3.climb.ac.uk/kraken2-microbial/hash.k2d
    wget https://refdb.s3.climb.ac.uk/kraken2-microbial/opts.k2d
    wget https://refdb.s3.climb.ac.uk/kraken2-microbial/taxo.k2d
    """
}
