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
include { DO_MLST                        } from '../subworkflows/local/do_mlst'

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

workflow UPDATE_CDC_PHOENIX_EXQC {
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

        CDIFF_CLADE(
            CREATE_INPUT_CHANNELS.out.combined_mlst, ASSET_CHECK.out.mlst_db
        )
        ch_versions = ch_versions.mix(CDIFF_CLADE.out.versions)

        CDIFF_PLASMID(
            CREATE_INPUT_CHANNELS.out.combined_mlst, ASSET_CHECK.out.mlst_db
        )
        ch_versions = ch_versions.mix(CDIFF_PLASMID.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        filtered_scaffolds_ch = CREATE_INPUT_CHANNELS.out.filtered_scaffolds.map{    meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}
        .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome      -> [[id:meta.id, project_id:meta.project_id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])

        CDIFF_GENES(
            filtered_scaffolds_ch, params.cdiff_db
        )
        ch_versions = ch_versions.mix(CDIFF_GENES.out.versions)

        CENTAR_CONSOLIDATER(
            CDIFF_GENES.out.gamma, CDIFF_CLADE.out.clade
        )
        ch_versions = ch_versions.mix(CENTAR_CONSOLIDATER.out.versions)

        /// need to figure out how to do this on a per directory 
        // Collecting the software versions
        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

    /*emit:
        mlst             = DO_MLST.out.checked_MLSTs
        amrfinder_output = AMRFINDERPLUS_RUN.out.report
        gamma_ar         = GAMMA_AR.out.gamma
        phx_summary      = GATHER_SUMMARY_LINES.out.summary_report*/
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