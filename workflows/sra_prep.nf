/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

/*
========================================================================================
    IMPORT SUBWORKFLOWS
========================================================================================
*/

//
// SUBWORKFLOW: Read in samplesheet/list, validate and stage input files
//

include { CREATE_SRA_SAMPLESHEET         } from '../modules/local/create_sra_samplesheet'
include { SRATOOLS_PREFETCH              } from '../modules/local/sratools/prefetch'
include { SRATOOLS_FASTERQDUMP           } from '../modules/local/sratools/fasterq'
include { ENTREZDIRECT_ESEARCH           } from '../modules/local/get_sra_metadata'
include { RENAME_SRA_FASTA               } from '../modules/local/rename_sra'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/

def add_seq_name(reads, metadata_csv) {
    // function is used to swap out SRR meta.id for the sample name from SRA
    def meta = [:] // create meta array
    meta.id = metadata_csv.readLines().get(1).split(',')[29] // This gives the metadata sample name from the SRA
    output_array = [ meta, reads]
    return output_array
}

def add_meta(folder) {
    def meta = [:] // create meta array
    meta.id = folder.toString().substring(folder.toString().lastIndexOf("/") + 1) // This gets the metadata sample name from the SRA, the +1 removes the /
    output_array = [ meta, folder]
    // Brief pause to keep the ncbi requests from hitting a limit
    sleep(5000) // sleep takes number as millisecond so this will pause for 5 seconds
    return output_array
}

file_count = 0
workflow SRA_PREP {
    take:
        sra_samplesheet //params.input_sra

    main:
        ch_versions = Channel.empty() // Used to collect the software versions
        outdir_path = Channel.fromPath(params.outdir, relative: true)

    SRATOOLS_PREFETCH (
        sra_samplesheet
    )
    ch_versions = ch_versions.mix(SRATOOLS_PREFETCH.out.versions)
    
    //flatten to have one folder go in at a time for parallel processing
    sra_folder_ch = SRATOOLS_PREFETCH.out.sra_folder.flatten().map{ it -> add_meta(it) }

    SRATOOLS_FASTERQDUMP (
        sra_folder_ch
    )
    ch_versions = ch_versions.mix(SRATOOLS_FASTERQDUMP.out.versions)

    //Getting the sample name from the metadata on SRR page
    ENTREZDIRECT_ESEARCH (
        sra_folder_ch, "sra"
    )
    ch_versions = ch_versions.mix(ENTREZDIRECT_ESEARCH.out.versions)

    // First, create channel with reads and metadata in it and make sure the right metadata goes with the right reads
    rename_sra_ch = SRATOOLS_FASTERQDUMP.out.reads.join(ENTREZDIRECT_ESEARCH.out.metadata_csv, by: [0])

    // Rename SRAs to have illumina style extension so PHX doesn't complain during input_check
    // Also, get sample name from metadata and rename fastq file
    RENAME_SRA_FASTA ( 
        // Create tuple with sample name in it for the tag in RENAME_SRA_FASTA
        rename_sra_ch.map{meta, reads, metadata_csv -> add_seq_name(reads, metadata_csv) }
    )
    ch_versions = ch_versions.mix(RENAME_SRA_FASTA.out.versions)

    // Gather fastqs and metadata to be sent to results folder
    CREATE_SRA_SAMPLESHEET (
        //removing the meta information because we are collecting everything so that part doesn't matter
        SRATOOLS_FASTERQDUMP.out.reads.map{meta, reads -> reads }.collect(), ENTREZDIRECT_ESEARCH.out.metadata_csv.map{meta, metadata_csv -> metadata_csv }.collect(), outdir_path
    )
    ch_versions = ch_versions.mix(CREATE_SRA_SAMPLESHEET.out.versions)

    emit:
        samplesheet = CREATE_SRA_SAMPLESHEET.out.csv  // stored in results folder
        versions    = ch_versions                    // channel: [ versions.yml ]
}
