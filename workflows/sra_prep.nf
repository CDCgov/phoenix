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
include { CONFIRM_DUPS                   } from '../subworkflows/local/confirm_dups'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/

def add_meta(folder) {
    def meta = [:] // create meta array
    meta.id = folder.toString().substring(folder.toString().lastIndexOf("/") + 1) - "_Folder" // This gets the metadata sample name from the SRA, the +1 removes the /
    def output_array = [ meta, folder]
    return output_array
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SRA_PREP {
    take:
        sra_samplesheet //params.input_sra

    main:
        ch_versions = Channel.empty() // Used to collect the software versions
        outdir_path = Channel.fromPath(params.outdir, relative: true)
        // slipt samplesheet to have one SRR come out at a time
        sra_id_ch = Channel.from(sra_samplesheet).splitCsv()

    // Fetch .sra data for each SRR number
    SRATOOLS_PREFETCH (
        sra_id_ch
    )
    ch_versions = ch_versions.mix(SRATOOLS_PREFETCH.out.versions)

    // Collect all output folders then flatten to have one folder go in at a time for parallel processing and metering of jobs so we don't hit NCBI's request limit
    sra_folder_ch = SRATOOLS_PREFETCH.out.sra_folder.map{ it -> add_meta(it) }

    // Get R1 and R2 files from .sra file
    SRATOOLS_FASTERQDUMP (
        sra_folder_ch
    )
    ch_versions = ch_versions.mix(SRATOOLS_FASTERQDUMP.out.versions)

    // Getting the metadata csv file for each SRR to get sample name
    ENTREZDIRECT_ESEARCH (
        sra_folder_ch
    )
    ch_versions = ch_versions.mix(ENTREZDIRECT_ESEARCH.out.versions)

    // First, create channel with reads and metadata in it and make sure the right metadata goes with the right reads
    combined_sra_ch = SRATOOLS_FASTERQDUMP.out.reads.join(ENTREZDIRECT_ESEARCH.out.metadata_csv, by: [0])

    // Figuring out if we should use SRR or sample names
    CONFIRM_DUPS (
        combined_sra_ch, ENTREZDIRECT_ESEARCH.out.metadata_csv
    )

    // Rename SRAs to have illumina style extension so PHX doesn't complain during input_check
    // Also, get sample name from metadata and rename fastq file
    RENAME_SRA_FASTA (
        CONFIRM_DUPS.out.rename_sra_out_ch
    )
    ch_versions = ch_versions.mix(RENAME_SRA_FASTA.out.versions)

    // Gather fastqs and metadata to be sent to results folder
    CREATE_SRA_SAMPLESHEET (
        RENAME_SRA_FASTA.out.renamed_reads.collect(), \
        //removing the meta information because we are collecting everything so that part doesn't matter
        ENTREZDIRECT_ESEARCH.out.metadata_csv.map{meta, metadata_csv -> metadata_csv }.collect(), \
        outdir_path, \
        params.use_sra
    )
    ch_versions = ch_versions.mix(CREATE_SRA_SAMPLESHEET.out.versions)

    emit:
        samplesheet = CREATE_SRA_SAMPLESHEET.out.csv  // stored in results folder
        versions    = ch_versions                    // channel: [ versions.yml ]
}

// Adding if/else for running on ICA
if (params.ica==false) {
    // do nothing, not running ICA and no erros occurred
} else if (params.ica==true) {
    workflow.onError { 
        // copy intermediate files + directories
        println("Getting intermediate files from ICA")
        ['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
        // return trace files
        println("Returning workflow run-metric reports from ICA")
        ['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'grep','"report"' ,'|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
    }
} else {
        error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods."
}