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

include { ASSET_CHECK                          } from '../modules/local/asset_check'
include { GAMMA as GAMMA_AR                    } from '../modules/local/gamma'
include { GET_TAXA_FOR_AMRFINDER               } from '../modules/local/get_taxa_for_amrfinder'
include { AMRFINDERPLUS_RUN                    } from '../modules/local/run_amrfinder'
include { CREATE_SUMMARY_LINE                  } from '../modules/local/phoenix_summary_line'
include { CREATE_AND_UPDATE_README             } from '../modules/local/updater/create_and_update_readme'
include { FETCH_FAILED_SUMMARIES               } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES                 } from '../modules/local/phoenix_summary'
include { GRIPHIN as GRIPHIN_NO_PUBLISH        } from '../modules/local/griphin'
include { GRIPHIN as GRIPHIN_NO_PUBLISH_CDC    } from '../modules/local/griphin'
include { UPDATE_GRIPHIN                       } from '../modules/local/updater/update_griphin'
include { UPDATE_GRIPHIN as UPDATE_CDC_GRIPHIN } from '../modules/local/updater/update_griphin'
include { SRST2_AR                             } from '../modules/local/srst2_ar'

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

def remove_last_element(input_ch) {
    return input_ch.remove(-1)
}

def add_meta(full_path_txt, griphin_ch) {
    def meta = [:] // create meta array
    def folder_path = full_path_txt.readLines().first()
    //meta.project_id = input_ch.getName().replaceAll("_GRiPHin.xlsx", "") // get file name without extention
    meta.project_id = folder_path
    def array = [ meta, griphin_ch ]
    return array
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow UPDATE_PHOENIX_WF {
    take:
        ch_input
        ch_input_indir
        ch_versions

    main:
        // Allow outdir to be relative
        outdir_path = Channel.fromPath(params.outdir, relative: true)

        CREATE_INPUT_CHANNELS (
            // False is to say if centar is to be included when creating input channels
            ch_input_indir, ch_input, false
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
            params.zipped_sketch, params.custom_mlstdb, []
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

        // combing fastp_trimd information with fairy check of reads to confirm there are reads after filtering
        trimd_reads_file_integrity_ch = CREATE_INPUT_CHANNELS.out.reads.join(CREATE_INPUT_CHANNELS.out.fairy_outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])
        .join(CREATE_INPUT_CHANNELS.out.phoenix_tsv_ch.map{meta, phoenix_tsv_ch -> [[id:meta.id, project_id:meta.project_id], meta.entry]}, by: [[0][0],[0][1]])

        // now we will split the channel into its true (busco present) and false (busco wasn't run with this dataset) elements
        busco_boolean_1ch = trimd_reads_file_integrity_ch.branch{ 
                    buscoTrue: it[3] == true
                    buscoFalse: it[3] == false}

        // Idenitifying AR genes in trimmed reads - using only datasets that were previously run with CDC_PHOENIX entry
        SRST2_AR (
            busco_boolean_1ch.buscoTrue.map{meta, reads, fairy_outcome, busco_boolean -> [meta, reads, fairy_outcome]}, "gene", params.ardb
        )
        ch_versions = ch_versions.mix(SRST2_AR.out.versions)

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

        // Extract the list of project folders from input_dir channel, combine with the actual path that is needed by GATHER_SUMMARY_LINES
        project_ids = CREATE_INPUT_CHANNELS.out.directory_ch.map { it[1] }.unique().collect()
            .map { dirs -> dirs.collectEntries { dir -> def dirName = dir.tokenize('/').last()
            return [dirName, dir]}}

        // Group the line_summaries by their project id and add in the full path for the project dir
        //the join is only to make sure SRST2 finished before going to the last step, we just combine and then kick it out
        summaries_ch = CREATE_SUMMARY_LINE.out.line_summary
            .join(SRST2_AR.out.gene_results.map{meta, gene_results -> [[id:meta.id, project_id:meta.project_id], gene_results]}, by: [[0][0],[0][1]]) 
            .map{meta, summaryline, gene_results -> [meta.project_id, summaryline]}.groupTuple(by: 0)
            .map { group -> 
                def (id, files) = group
                id2 = id.split('/')[-1] // Remove full path and cut just the project id for use
                def dirPath = project_ids.value[id2]  // Retrieve the matching directory path
                return [id, dirPath, files]}

        //println("OMG")
        //summaries_ch.view()

        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            summaries_ch.map{ project_id, dir, summary_lines -> summary_lines}, summaries_ch.map{ project_id, dir, summary_lines -> dir}, true
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        // check if the incoming file was created with CDC_PHOENIX or PHOENIX
        buscoExists_ch = summaries_ch.map{ project_id, dir, summary_lines -> 
                def busco_boolean = summary_lines.head().text.contains('BUSCO')
                return [ summary_lines, dir, busco_boolean ]}

        // now we will split the channel into its true and false elements for if we need to run the CDC version of GRIPHIN or not. 
        busco_boolean_ch = buscoExists_ch.branch{
                    buscoTrue: it[2] == true 
                    buscoFalse: it[2] == false}

        // for samples that were created with entry CDC_PHX
        GRIPHIN_NO_PUBLISH_CDC (
            busco_boolean_ch.buscoTrue.map{ summary_line, dir, busco_boolean -> summary_line}, \
            CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, \
            busco_boolean_ch.buscoTrue.map{ summary_line, dir, busco_boolean -> dir.toString()}, params.coverage, false, false, true, false, false
        )
        ch_versions = ch_versions.mix(GRIPHIN_NO_PUBLISH_CDC.out.versions)

        // for samples that were created with entry PHX
        GRIPHIN_NO_PUBLISH (
            busco_boolean_ch.buscoFalse.map{ summary_line, dir, busco_boolean -> summary_line}, \
            CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, \
            busco_boolean_ch.buscoFalse.map{ summary_line, dir, busco_boolean -> dir.toString()}, params.coverage, true, false, true, false, false
        )
        ch_versions = ch_versions.mix(GRIPHIN_NO_PUBLISH.out.versions)

        // bring all the griphins into one channel and pass one at a time to the UPDATE_GRIPHIN process
<<<<<<< HEAD
        griphin_reports_ch = GRIPHIN_NO_PUBLISH.out.griphin_report.collect().ifEmpty([]).combine(GRIPHIN_NO_PUBLISH_CDC.out.griphin_report.collect().ifEmpty([])).flatten()
        //CREATE_INPUT_CHANNELS.out.griphin_excel_ch.view()
        //griphin_reports_ch.view()
        //CREATE_INPUT_CHANNELS.out.directory_ch.view()

        //map{it -> crush_meta(it)}

        // join old and new griphins for combining
        griphins_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch.map{it -> crush_meta(it)}.map{   meta, griphin_excel_ch -> [[project_id:meta.project_id], griphin_excel_ch]}\
        .join(griphin_reports_ch.map{it -> add_meta(it)}.map{                                     meta, griphin_report   -> [[project_id:meta.project_id], griphin_report]}, by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{                                         meta, directory_ch     -> [[project_id:meta.project_id], directory_ch]},   by: [[0][0],[0][1]])
=======
        griphin_reports_ch = GRIPHIN_NO_PUBLISH.out.griphin_report.collect().ifEmpty([]).combine(GRIPHIN_NO_PUBLISH_CDC.out.griphin_report.collect().ifEmpty([])).flatten().collate(2)\
                                .map{path_txt, griphin_report -> add_meta(path_txt, griphin_report)}

        // join old and new griphins for combining
        griphins_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch.map{meta, griphin_excel_ch -> [[project_id:meta.project_id], griphin_excel_ch]}\
        .join(griphin_reports_ch.map{meta, griphin_report   -> [[project_id:meta.project_id], griphin_report]}, by: [[0][0],[0][1]])\
        .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{            meta, directory_ch     -> [[project_id:meta.project_id], directory_ch]},   by: [[0][0],[0][1]])
>>>>>>> f579344 (updater fix)
        //.join(CREATE_INPUT_CHANNELS.out.valid_samplesheet.map {      meta, valid_samplesheet      -> [[project_id:meta.project_id], valid_samplesheet]}, by: [[0][0],[0][1]])

        //griphins_ch.view()

        // combine griphin files, the new one just created and the old one that was found in the project dir. 
        UPDATE_GRIPHIN (
            griphins_ch.map{ meta, old_excel, new_excel, directory -> [ old_excel, new_excel ] }, 
            griphins_ch.map{ meta, old_excel, new_excel, directory -> directory },
            griphins_ch.map{ meta, old_excel, new_excel, directory -> meta.project_id.toString().split('/')[-1] },
            //griphins_ch.map{ meta, excel, report, directory, samplesheet -> samplesheet },
            [],
            params.coverage
        )
        ch_versions = ch_versions.mix(UPDATE_GRIPHIN.out.versions)

        /// need to figure out how to do this on a per directory 
        // Collecting the software versions
        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

    emit:
        mlst             = DO_MLST.out.checked_MLSTs
        amrfinder_output = AMRFINDERPLUS_RUN.out.report
        gamma_ar         = GAMMA_AR.out.gamma
        phx_summary      = GATHER_SUMMARY_LINES.out.summary_report
        //output for phylophoenix
        griphin_tsv      = UPDATE_GRIPHIN.out.griphin_tsv_report
        griphin_excel    = UPDATE_GRIPHIN.out.griphin_report
        //dir_samplesheet  = UPDATE_CDC_GRIPHIN.out.converted_samplesheet
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
