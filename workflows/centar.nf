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

// Info required for completion email and summary
def multiqc_report = []

/*
========================================================================================
    IMPORT LOCAL MODULES
========================================================================================
*/

include { ASSET_CHECK                                         } from '../modules/local/asset_check'
include { CREATE_SUMMARY_LINE                                 } from '../modules/local/phoenix_summary_line'
include { FETCH_FAILED_SUMMARIES                              } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES as CENTAR_GATHER_SUMMARY_LINES } from '../modules/local/phoenix_summary' // calling it centar so output can stay together. 
include { GRIPHIN as CENTAR_GRIPHIN                           } from '../modules/local/griphin' // calling it centar so output can stay together.
include { UPDATE_GRIPHIN                                      } from '../modules/local/updater/update_griphin' 

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/
include { CREATE_INPUT_CHANNELS          } from '../subworkflows/local/create_input_channels'
include { GENERATE_PIPELINE_STATS_WF     } from '../subworkflows/local/generate_pipeline_stats'
include { CENTAR_SUBWORKFLOW             } from '../subworkflows/local/centar_steps'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { CUSTOM_DUMPSOFTWAREVERSIONS as CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main_centar'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/


def filter_out_meta(input_ch) {
    // Create an empty list to store the filtered items
    def filteredList = []
    def indir = input_ch[-1]
    // Iterate over the input list with indices
    input_ch.eachWithIndex { item, index ->
        // Add items at even indices to the filtered list
        if (index % 2 == 1) {
            filteredList << item
        }
    }
    return [ filteredList, indir ]
}

def create_meta(input_ch) {
    def meta = [:] // create meta array
    meta.project_id = input_ch.getName().replaceAll("_GRiPHin_Summary.xlsx", "")
    array = [ meta, input_ch ]
    return array
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow RUN_CENTAR {
    take:
        ch_input
        ch_input_indir
        ch_versions
        outdir_path

    main:
        CREATE_INPUT_CHANNELS (
            ch_input_indir, ch_input, true
        )
        ch_versions = ch_versions.mix(CREATE_INPUT_CHANNELS.out.versions)

        //unzip any zipped databases
        ASSET_CHECK (
            params.zipped_sketch, params.custom_mlstdb, []
        )
        ch_versions = ch_versions.mix(ASSET_CHECK.out.versions)

        // run centar workflow
        CENTAR_SUBWORKFLOW (
            CREATE_INPUT_CHANNELS.out.combined_mlst,
            CREATE_INPUT_CHANNELS.out.fairy_outcome,
            CREATE_INPUT_CHANNELS.out.filtered_scaffolds,
            ASSET_CHECK.out.mlst_db,
            CREATE_INPUT_CHANNELS.out.taxonomy
        )
        ch_versions = ch_versions.mix(CENTAR_SUBWORKFLOW.out.versions)

        // get summary lines and directory information to make sure all samples for a particular project folder stay together. 
        summaries_ch = CREATE_INPUT_CHANNELS.out.line_summary.map{meta, line_summary -> [[project_id:meta.project_id], line_summary] }
            .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{     meta, dir          -> [[project_id:meta.project_id], dir]}, by: [0])
            .map { it -> filter_out_meta(it) }

        // Combining sample summaries into final report
        CENTAR_GATHER_SUMMARY_LINES (
            summaries_ch.map{ summary_line, dir -> summary_line}, summaries_ch.map{ summary_line, dir -> dir}, true
        )
        ch_versions = ch_versions.mix(CENTAR_GATHER_SUMMARY_LINES.out.versions)

        // collect centar output and summary lines to make the griphin report
        griphin_input = CREATE_INPUT_CHANNELS.out.line_summary.map{meta, line_summary      -> [[project_id:meta.project_id], line_summary] }.groupTuple(by: [0])\
            .join(CENTAR_SUBWORKFLOW.out.consolidated_centar.map{  meta, consolidated_file -> [[project_id:meta.project_id], consolidated_file]}.groupTuple(by: [0]), by: [0])\
            .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{      meta, dir               -> [[project_id:meta.project_id], dir]}, by: [0])

       // separate the summary and centar files from dir
        griphin_input_ch = griphin_input.map{meta, summary_line, centar_files, dir -> [summary_line, centar_files].flatten()}
        // get project dir - indir
        griphin_dir_path = griphin_input.map{meta, summary_line, centar_files, dir -> [dir]}

        //create GRiPHin report
        CENTAR_GRIPHIN (
            griphin_input_ch, CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, griphin_dir_path, params.coverage, false, false, false, false, CENTAR_SUBWORKFLOW.out.run_centar_in_griphin
        )
        ch_versions = ch_versions.mix(CENTAR_GRIPHIN.out.versions)

        if (params.combine_griphins == false) { // don't run if combined workflow is going to run
            /*if (params.outdir != "${launchDir}/phx_output"){ // when no outdir is passed
                // to be able to create software_versions.yml 
                software_versions_ch = ch_versions.unique().collectFile(name: 'collated_versions.yml') // combine with CENTAR_GRIPHIN.out to ensure this runs after
                .combine(params.outdir).map{version, outdir -> [outdir.toString(), version]}
            } else {*/
                // to be able to create software_versions.yml 
                software_versions_ch = ch_versions.unique().collectFile(name: 'collated_versions.yml') // combine with CENTAR_GRIPHIN.out to ensure this runs after
                .combine(CENTAR_GRIPHIN.out.griphin_report).map{version, meta_file, griphin -> [meta_file.splitText().first().toString().trim(), version]}
            //}

            /////////////////////////////// need to figure out how to do this on a per directory 
            // Collecting the software versions
            CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS (
                software_versions_ch.map{meta_file, version -> version},
                software_versions_ch.map{projectDir, version -> projectDir}
            )
        } else {
            software_versions_ch = ch_versions
        }

    emit:
        //output for phylophoenix
        griphins_tsv    = CENTAR_GRIPHIN.out.griphin_tsv_report
        griphins_excel  = CENTAR_GRIPHIN.out.griphin_report
        dir_samplesheet = CENTAR_GRIPHIN.out.converted_samplesheet
        valid_samplesheet = CREATE_INPUT_CHANNELS.out.valid_samplesheet
        ch_versions     = software_versions_ch
}

if (params.combine_griphins == false) { // don't run if combined workflow is going to run
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
                    NfcoreTemplate.centar_email(workflow, params, summary_params, projectDir, log, multiqc_report)
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
}

/*
========================================================================================
    THE END
========================================================================================
*/