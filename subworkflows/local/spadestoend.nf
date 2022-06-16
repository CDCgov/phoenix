include { GAMMA } from '../modules/nf-core/modules/gamma/main'
include { FASTP } from '../modules/nf-core/modules/fastp/main'
include { SPADES } from '../modules/nf-core/modules/spades/main'
include { QUAST } from '../modules/nf-core/modules/quast/main'
include { FASTANI } from '../modules/nf-core/modules/fastani/main'
include { MASH_DIST } from '../modules/nf-core/modules/mash/dist/main'
include { MLST } from '../modules/nf-core/modules/mlst/main'
include { KRAKEN2 } from '../modules/nf-core/modules/kraken2/main'
include { PROKKA } from '../modules/nf-core/modules/prokka/main'
include { BUSCO } from '../modules/nf-core/modules/busco/main'
include { KRONA } from '../modules/nf-core/modules/krona/main'

//local module
include { KRAKEN2_DB } from '../modules/local/kraken2db'

// database parameter checks

if(params.gamma_db){
    Channel
        .fromPath( "${params.gamma_db}" )
        .set { ch_gamma }
} else {
    ch_gamma = Channel.empty()
}


workflow SPADES_TO_END {
    
    KRONA_KRONADB ( ) 

    KRONA_KTIMPORTTAXONOMY ( Need map, KRONA_KRONADB.out.db, taxes.csv?  ) 

    SPADES ( FASTP.out.reads, directry/file for aa HMMS for guided mode?)

    QUAST( SPADES.out.scaffolds, SPADES.out.contigs, true, SPADES.out.gfa, true )

    FASTANI ( SPADES.out.scaffolds, SPADES.out.scaffolds , reference file for query)

    MASHTREE ( SPADES.out.scaffolds )

    MLST ( SPADES.out.scaffolds ) 

    GAMMA ( SPADES.out.scaffolds, ch_gamma )

    KRAKEN2_DB (SPADES.out.scaffolds, KRAKEN2_DB.out.db )

    KRAKEN2 (SPADES.out.scaffolds, KRAKEN2_DB.out.db) 

    PROKKA ( SPADES.out.scaffolds )

    BUSCO () 
}