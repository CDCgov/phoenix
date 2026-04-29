/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

/*
========================================================================================
    SETUP
========================================================================================
*/

/*
========================================================================================
    IMPORT LOCAL MODULES
========================================================================================
*/
include { FILE_RENAME                               } from '../modules/local/file_rename'
include { UPDATE_GRIPHIN                            } from '../modules/local/updater/update_griphin'
include { UPDATE_GRIPHIN as UPDATE_GRIPHIN_ONLY_TWO } from '../modules/local/updater/update_griphin'

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK                    } from '../subworkflows/local/input_check'
include { CREATE_INPUT_CHANNELS          } from '../subworkflows/local/create_input_channels'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS as COMBINE_CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { META_CUSTOM_DUMPSOFTWAREVERSIONS                                    } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main_species_specific'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/

def add_meta(input_ch) {
    def meta = [:] // create meta array
    meta.project_id = input_ch.getName().replaceAll("_GRiPHin.xlsx", "") // get file name without extention
    def array = [ meta, input_ch ]
    return array
}

def findGRiPHinSummaryFiles(path, extension) {
        def files = file(path[0]).listFiles() // List all files in the directory
        files.findAll { it.name.endsWith('_GRiPHin_Summary.'+ extension) } // Find all matching files
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow COMBINE_GRIPHINS_WF {
    take:
        ch_input
        ch_versions

    main:

        // SUBWORKFLOW: Read in samplesheet/list, validate and stage input files
        INPUT_CHECK (
            ch_input, true
        )
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        // Rename GRiPHin files so there isn't collision when combining
        FILE_RENAME (
            INPUT_CHECK.out.griphins
        )
        ch_versions = ch_versions.mix(FILE_RENAME.out.versions)

        def griph_out_path
        // if there is no custom output directory for griphins, use the main outdir
        if (params.outdir == "${launchDir}/phx_output" && params.griphin_out != "${launchDir}") {
            griph_out_path = params.griphin_out
        } else {
            // Allow outdir to be relative
            griph_out_path = Channel.fromPath(params.outdir, relative: true)
        }

        griphins_ch = FILE_RENAME.out.renamed_griphins.toList()combine(griph_out_path).map{ files, outdir_path ->
                def meta = [:]
                meta.project_id = outdir_path.toString().split('/')[-1]
                meta.full_project_id = outdir_path
                return [meta, files, files.size()]
            }.branch{
                exactly_two: it[2] == 2
                more_than_two: it[2] > 2
            }

        exactly_two_ch = griphins_ch.exactly_two.map{ meta, fileList, size -> [ meta, [fileList[0], fileList[1]] ] }
        more_than_two_ch = griphins_ch.more_than_two.map{ meta, fileList, size -> [ meta, fileList ] }

        // combine griphin files, the new one just created and the old one that was found in the project dir. 
        UPDATE_GRIPHIN_ONLY_TWO (
            exactly_two_ch.map{ meta, old_excel, new_excel -> [old_excel, new_excel] },
            exactly_two_ch.map{ meta, old_excel, new_excel -> meta.full_project_id },
            INPUT_CHECK.out.valid_samplesheet,
            params.coverage,
            params.bldb,
            true,
            exactly_two_ch.map{ meta, old_excel, new_excel -> meta.project_id }
        )
        ch_versions = ch_versions.mix(UPDATE_GRIPHIN_ONLY_TWO.out.versions)

        // combine griphin files, the new one just created and the old one that was found in the project dir. 
        UPDATE_GRIPHIN (
            more_than_two_ch.map{ meta, excels -> excels },
            more_than_two_ch.map{ meta, excels -> meta.full_project_id },
            INPUT_CHECK.out.valid_samplesheet,
            params.coverage,
            params.bldb,
            true,
            // Some boolean to indicate if this is a species specific entry point
            more_than_two_ch.map{ meta, excels -> meta.project_id }
        )
        ch_versions = ch_versions.mix(UPDATE_GRIPHIN.out.versions)

        griphin_tsv_report = UPDATE_GRIPHIN.out.griphin_tsv_report.ifEmpty([]).mix(UPDATE_GRIPHIN_ONLY_TWO.out.griphin_tsv_report.ifEmpty([]))
        griphin_report = UPDATE_GRIPHIN.out.griphin_report.ifEmpty([]).mix(UPDATE_GRIPHIN_ONLY_TWO.out.griphin_report.ifEmpty([]))

        /// need to figure out how to do this on a per directory
        // Collecting the software versions
        COMBINE_CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

    emit:
        //output for phylophoenix
        griphin_tsv      = griphin_tsv_report
        griphin_excel    = griphin_report
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

// Adding if/else for running on ICA
if (params.ica==false) {
    workflow.onComplete {
        if (count == 0) {
            if (params.email || params.email_on_fail) {
                NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
            }
            NfcoreTemplate.summary(workflow, params, log)
            count++
        }
    }
} else if (params.ica==true) {
    workflow.onComplete { 
        if (workflow.success) {
            println("Pipeline Completed Successfully")
            NfcoreTemplate.summary(workflow, params, log)
            System.exit(0)
        } else {
            println("Pipeline Completed with Errors")
            NfcoreTemplate.summary(workflow, params, log)
            System.exit(1)
        }
    }
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

/*
========================================================================================
    THE END
========================================================================================
*/
