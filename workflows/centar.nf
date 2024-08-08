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

include { ASSET_CHECK                    } from '../modules/local/asset_check'
include { CDIFF_CLADE                    } from '../modules/local/centar/cdiff_clade'
include { CDIFF_PLASMID                  } from '../modules/local/centar/cdiff_plasmid'
include { GAMMA as CDIFF_GENES           } from '../modules/local/gamma'
include { CENTAR_CONSOLIDATER            } from '../modules/local/centar/centar_consolidater'
include { CREATE_SUMMARY_LINE            } from '../modules/local/phoenix_summary_line'
include { FETCH_FAILED_SUMMARIES         } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES           } from '../modules/local/phoenix_summary'
include { GRIPHIN as GRIPHIN_COPY        } from '../modules/local/griphin'
include { UPDATE_GRIPHIN                 } from '../modules/local/updater/update_griphin'

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

include { CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

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

    main:
        // Allow outdir to be relative
        outdir_path = Channel.fromPath(params.outdir, relative: true)
        // Allow relative paths for krakendb argument
        kraken2_db_path  = Channel.fromPath(params.kraken2db, relative: true)

        CREATE_INPUT_CHANNELS (
            ch_input_indir, ch_input, ch_versions
        )
        ch_versions = ch_versions.mix(CREATE_INPUT_CHANNELS.out.versions)

        //unzip any zipped databases
        ASSET_CHECK (
            params.zipped_sketch, params.custom_mlstdb, kraken2_db_path
        )
        ch_versions = ch_versions.mix(ASSET_CHECK.out.versions)

        CENTAR_SUBWORKFLOW (
            CREATE_INPUT_CHANNELS.out.combined_mlst,
            CREATE_INPUT_CHANNELS.out.fairy_outcome,
            CREATE_INPUT_CHANNELS.out.filtered_scaffolds,
            ASSET_CHECK.out.mlst_db
        )
        ch_versions = ch_versions.mix(CENTAR_SUBWORKFLOW.out.versions)

        // Extract the list of project folders from input_dir channel
        project_ids = CREATE_INPUT_CHANNELS.out.directory_ch.map { it[1] }.collect().toList()

        // Collect all the summary files and filter to keep those that have the correct project_id
        summaries_ch = CREATE_SUMMARY_LINE.out.line_summary.collect().filter { line_summary ->
            def project_id = line_summary[0].project_id
            project_ids.contains(project_id)}
            .join(CREATE_INPUT_CHANNELS.out.directory_ch, by: [[0][0],[0][1]])
            .map { it -> filter_out_meta(it) }// Convert to a list to use indexing

        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            summaries_ch.map{ summary_line, dir -> summary_line}, summaries_ch.map{ summary_line, dir -> dir}, true
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        GRIPHIN_COPY (
            summaries_ch.map{ summary_line, dir -> summary_line}, \
            CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, \
            summaries_ch.map{ summary_line, dir -> dir.toString()}, params.coverage, true, false
        )
        ch_versions = ch_versions.mix(GRIPHIN_COPY.out.versions)

        // combine by project_id - > this will need to be adjusted in CDC_PHOENIX and PHOENIX
        griphins_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch.map{         meta, griphin_excel_ch -> [[project_id:meta.project_id], griphin_excel_ch]}\
        .join(GRIPHIN_COPY.out.griphin_report.map{ it -> create_meta(it)}.map{meta, griphin_report   -> [[project_id:meta.project_id], griphin_report]}, by: [0])\
        .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{                     meta, directory_ch     -> [[project_id:meta.project_id], directory_ch]},   by: [0])

        griphins_ch.view()

        UPDATE_GRIPHIN (
            griphins_ch, params.coverage
        )
        ch_versions = ch_versions.mix(UPDATE_GRIPHIN.out.versions)

        /// need to figure out how to do this on a per directory 
        // Collecting the software versions
        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

        //
        // MODULE: MultiQC
        //
        workflow_summary    = WorkflowPhoenix.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    emit:
        amrfinder_output = CENTAR_CONSOLIDATER.out.report*/

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