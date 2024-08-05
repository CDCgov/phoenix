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
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES
========================================================================================
*/

include { ASSET_CHECK                    } from '../modules/local/asset_check'
include { GAMMA as GAMMA_AR              } from '../modules/local/gamma'
include { GET_TAXA_FOR_AMRFINDER         } from '../modules/local/get_taxa_for_amrfinder'
include { AMRFINDERPLUS_RUN              } from '../modules/local/run_amrfinder'
include { CREATE_SUMMARY_LINE            } from '../modules/local/phoenix_summary_line'
include { CREATE_AND_UPDATE_README       } from '../modules/local/updater/create_and_update_readme'
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

include { MULTIQC                      } from '../modules/nf-core/modules/multiqc/main'
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

        // combine directory and pipeline_info - make sure that id and project folder match
        files_to_update_ch = CREATE_INPUT_CHANNELS.out.directory_ch.join(CREATE_INPUT_CHANNELS.out.pipeline_info, by: [[0][0],[0][1]])

        CREATE_AND_UPDATE_README (
            files_to_update_ch,
            workflow.manifest.version,
            params.custom_mlstdb,
            params.ardb,
            params.amrfinder_db
        )
        ch_versions = ch_versions.mix(CREATE_AND_UPDATE_README.out.versions)

        //unzip any zipped databases
        ASSET_CHECK (
            params.zipped_sketch, params.custom_mlstdb, kraken2_db_path
        )
        ch_versions = ch_versions.mix(ASSET_CHECK.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        filtered_scaffolds_ch = CREATE_INPUT_CHANNELS.out.filtered_scaffolds.map{    meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}
        .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome      -> [[id:meta.id, project_id:meta.project_id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])

        // Running gamma to identify AR genes in scaffolds
        GAMMA_AR (
            filtered_scaffolds_ch, params.ardb
        )
        ch_versions = ch_versions.mix(GAMMA_AR.out.versions)

        // Perform MLST steps on isolates (with srst2 on internal samples)
        // 3rd input would normally be FASTP_TRIMD.out.reads, but since we have no reads we just pass something to get the meta_id and create an empty channel
        // This is necessary as we need the meta.id information for the mid_srst2_ch in the subworkflow.
        DO_MLST (
            CREATE_INPUT_CHANNELS.out.filtered_scaffolds, \
            CREATE_INPUT_CHANNELS.out.fairy_outcome, \
            CREATE_INPUT_CHANNELS.out.reads, \
            CREATE_INPUT_CHANNELS.out.taxonomy, \
            ASSET_CHECK.out.mlst_db, \
            true, "update"
        )
        ch_versions = ch_versions.mix(DO_MLST.out.versions)

        // Create file that has the organism name to pass to AMRFinder
        GET_TAXA_FOR_AMRFINDER (
            CREATE_INPUT_CHANNELS.out.taxonomy, false
        )
        ch_versions = ch_versions.mix(GET_TAXA_FOR_AMRFINDER.out.versions)

        // Combining taxa and scaffolds to run amrfinder and get the point mutations.
        amr_channel = CREATE_INPUT_CHANNELS.out.filtered_scaffolds.map{          meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}\
        .join(GET_TAXA_FOR_AMRFINDER.out.amrfinder_taxa.splitCsv(strip:true).map{meta, amrfinder_taxa     -> [[id:meta.id, project_id:meta.project_id], amrfinder_taxa ]}, by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.prokka_faa.map{                          meta, prokka_faa         -> [[id:meta.id, project_id:meta.project_id], prokka_faa ]},     by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.prokka_gff.map{                          meta, prokka_gff         -> [[id:meta.id, project_id:meta.project_id], prokka_gff ]},     by: [[0][0],[0][1]])

        // Run AMRFinder
        AMRFINDERPLUS_RUN (
            amr_channel, params.amrfinder_db
        )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

         // Combining output based on meta.id to create summary by sample -- is this verbose, ugly and annoying? yes, if anyone has a slicker way to do this we welcome the input.
        line_summary_ch = CREATE_INPUT_CHANNELS.out.fastp_total_qc.map{meta, fastp_total_qc  -> [[id:meta.id, project_id:meta.project_id], fastp_total_qc]}\
        .join(DO_MLST.out.checked_MLSTs.map{                           meta, checked_MLSTs   -> [[id:meta.id, project_id:meta.project_id], checked_MLSTs]},  by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.gamma_hv.map{                  meta, gamma_hv        -> [[id:meta.id, project_id:meta.project_id], gamma_hv]},       by: [[0][0],[0][1]])\
        .join(GAMMA_AR.out.gamma.map{                                  meta, gamma           -> [[id:meta.id, project_id:meta.project_id], gamma]},          by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.gamma_pf.map{                  meta, gamma_pf        -> [[id:meta.id, project_id:meta.project_id], gamma_pf]},       by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.quast_report.map{              meta, quast_report    -> [[id:meta.id, project_id:meta.project_id], quast_report]},   by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.assembly_ratio.map{            meta, assembly_ratio  -> [[id:meta.id, project_id:meta.project_id], assembly_ratio]}, by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.synopsis.map{                  meta, synopsis        -> [[id:meta.id, project_id:meta.project_id], synopsis]},       by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.taxonomy.map{                  meta, taxonomy        -> [[id:meta.id, project_id:meta.project_id], taxonomy]},       by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.k2_bh_summary.map{             meta, k2_bh_summary   -> [[id:meta.id, project_id:meta.project_id], k2_bh_summary]},  by: [[0][0],[0][1]])\
        .join(AMRFINDERPLUS_RUN.out.report.map{                        meta, report          -> [[id:meta.id, project_id:meta.project_id], report]},         by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.ani_best_hit.map{              meta, ani_best_hit    -> [[id:meta.id, project_id:meta.project_id], ani_best_hit]},   by: [[0][0],[0][1]])

        // Generate summary per sample
        CREATE_SUMMARY_LINE (
            line_summary_ch
        )
        ch_versions = ch_versions.mix(CREATE_SUMMARY_LINE.out.versions)

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
        mlst             = DO_MLST.out.checked_MLSTs
        amrfinder_output = AMRFINDERPLUS_RUN.out.report
        gamma_ar         = GAMMA_AR.out.gamma
        phx_summary      = GATHER_SUMMARY_LINES.out.summary_report
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
