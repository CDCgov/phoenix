//
// SUBWORKFLOW: Read in samplesheet/list, validate and stage input files
//

include { SRA_SAMPLESHEET_CHECK          } from '../../modules/local/sra_samplesheet_check'
include { SRATOOLS_PREFETCH              } from '../../modules/local/sratools/prefetch'
include { SRATOOLS_FASTERQDUMP           } from '../../modules/local/sratools/fasterq'
include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { RENAME_SRA_FASTA               } from '../../modules/local/rename_sra'

workflow GET_SRA {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:

    ch_versions     = Channel.empty() // Used to collect the software versions

    SRATOOLS_PREFETCH (
        params.new_samplesheet
    )
    ch_versions = ch_versions.mix(SRATOOLS_PREFETCH.out.versions)

    SRATOOLS_FASTERQDUMP (
        SRATOOLS_PREFETCH.out.sra
    )
    ch_versions = ch_versions.mix(SRATOOLS_FASTERQDUMP.out.versions)

    //force delay parallelization to allow fastqs to be sent to results folder
    SRA_SAMPLESHEET_CHECK (
        SRATOOLS_PREFETCH.out.samplesheet, SRATOOLS_FASTERQDUMP.out.versions
    )
    ch_versions = ch_versions.mix(SRA_SAMPLESHEET_CHECK.out.versions)

    //Rename SRAs with input declared to stop immediate parallel processing
    RENAME_SRA_FASTA ( 
        SRATOOLS_FASTERQDUMP.out.versions
    )

    // Collecting the software versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    samplesheet        = SRA_SAMPLESHEET_CHECK.out.csv

    emit:
    samplesheet                       // stored in results folder
    versions = ch_versions           // channel: [ versions.yml ]
    
}