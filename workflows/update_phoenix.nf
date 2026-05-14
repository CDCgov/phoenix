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
include { SHIGAPASS                            } from '../modules/local/shigapass'
include { CHECK_SHIGAPASS_TAXA                 } from '../modules/local/check_shigapass_taxa'
include { CALCULATE_ASSEMBLY_RATIO             } from '../modules/local/assembly_ratio'
include { GAMMA as GAMMA_AR                    } from '../modules/local/gamma'
include { GAMMA as GAMMA_HV                    } from '../modules/local/gamma'
include { GAMMA_S as GAMMA_PF                  } from '../modules/local/gammas'
include { GET_TAXA_FOR_AMRFINDER               } from '../modules/local/get_taxa_for_amrfinder'
include { AMRFINDERPLUS_RUN                    } from '../modules/local/run_amrfinder'
include { CREATE_SUMMARY_LINE                  } from '../modules/local/phoenix_summary_line'
include { CREATE_AND_UPDATE_README             } from '../modules/local/updater/create_and_update_readme'
include { FETCH_FAILED_SUMMARIES               } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES                 } from '../modules/local/phoenix_summary'
include { GRIPHIN as GRIPHIN_NO_PUBLISH        } from '../modules/local/griphin'
include { GRIPHIN as GRIPHIN_PUBLISH           } from '../modules/local/griphin'
include { UPDATE_GRIPHIN                       } from '../modules/local/updater/update_griphin'
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

include { META_CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main_species_specific'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/
def remove_last_element(input_ch) {
    return input_ch.remove(-1)
}

def add_meta(full_path_txt, griphin_ch) {
    def meta = [:] // create meta array
    def folder_path = full_path_txt.readLines().first()
    //meta.project_id = input_ch.getName().replaceAll("_GRiPHin.xlsx", "") // get file name without extention
    meta.project_id = folder_path.toString().split('/')[-1]
    def array = [ meta, griphin_ch ]
    return array
}

def get_only_taxa(input_ch){ 
    def genus = ""
    def species = ""
    input_ch[1].eachLine { line ->
        if (line.startsWith("G:")) {
            def parts = line.split(":")[1].trim().split('\t')
            genus = parts.size() > 1 ? parts[1] : parts[0]
        } else if (line.startsWith("s:")) {
            def parts = line.split(":")[1].trim().split('\t')
            species = parts.size() > 1 ? parts[1] : parts[0]
        }
    }
    return [ "$genus" ]
}

def get_taxa(input_ch){ 
        def genus = ""
        def species = ""
        input_ch[1].eachLine { line ->
            if (line.startsWith("G:")) {
                def parts = line.split(":")[1].trim().split('\t')
                genus = parts.size() > 1 ? parts[1] : parts[0]
            } else if (line.startsWith("s:")) {
                def parts = line.split(":")[1].trim().split('\t')
                species = parts.size() > 1 ? parts[1] : parts[0]
            }
        }
        //return ["$genus $species", input_ch[0], input_ch[1]]
        return ["$genus", input_ch[0], input_ch[1]]
}

def pad_scaffold_samples = { source_ch, scaffold_anchor_ch, scaffold_pad_ch, empty_ch ->
    source_ch
        .join(scaffold_anchor_ch, by: [0])
        .map { meta, value, scaffolds -> [meta, value] }
        .ifEmpty(empty_ch)
        .mix(scaffold_pad_ch)
}

// Merge two channels and keep unique per sample
def merge_unique = { ch1, ch2 ->
    ch1
        .concat(ch2)
        .unique { meta, file -> [meta.id, meta.project_id] }
}

def filter_to_update_ids(source_ch, update_ids_ch) {
    return source_ch
        .combine(update_ids_ch)
        .filter{ meta, value, update_ids -> update_ids.contains(meta.id) }
        .map{ meta, value, update_ids -> [meta, value] }
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow UPDATE_PHOENIX_WF {
    take:
        ch_samplesheet_input
        ch_input_indir
        ch_versions

    main:
        // Allow outdir to be relative
        def outdir_path = Channel.fromPath(params.outdir, relative: true, type: 'dir')
        def empty_ch = Channel.empty()  // resolve ONCE, reuse everywhere
        def wf_version = workflow.manifest.version

        CREATE_INPUT_CHANNELS (
            // False is to say if centar is to be included when creating input channels
            ch_input_indir, ch_samplesheet_input, false
        )
        ch_versions = ch_versions.mix(CREATE_INPUT_CHANNELS.out.versions)

        CREATE_INPUT_CHANNELS.out.update_pipeline_info_isolate
            .collect()

        //unzip any zipped databases
        ASSET_CHECK (
            params.zipped_sketch, params.custom_mlstdb, [], params.clia_amrfinder_db
        )
        ch_versions = ch_versions.mix(ASSET_CHECK.out.versions)

        // Build sentinel ID list — gamma_ar is already date-filtered so it only
        // contains isolates that genuinely need updating.
        def orange = '\033[38;5;208m'
        def blue   = '\033[38;5;39m'
        def reset  = '\033[0m'

        isolates_needing_update = CREATE_INPUT_CHANNELS.out.sample_needs_update_ch
            .filter { meta, flag -> flag == true }
            .map { meta, flag -> meta.id }
            .collect()
            .map { ids -> [ids] }  // wrap to keep as single list element during combine
            .ifEmpty([["__EMPTY__"]])

        scaffolds_for_update_ch = CREATE_INPUT_CHANNELS.out.filtered_scaffolds
            .combine(isolates_needing_update)
            .filter{ meta, value, update_ids -> update_ids.contains(meta.id) }
            .map{ meta, value, update_ids -> [meta, value] }

        reads_for_update_ch = CREATE_INPUT_CHANNELS.out.reads
            .combine(isolates_needing_update)
            .filter{ meta, reads, update_ids -> update_ids.contains(meta.id) }
            .map{ meta, reads, update_ids -> [meta, reads] }

        CREATE_INPUT_CHANNELS.out.filtered_scaffolds
            .map{ meta, scaffolds -> meta.id }
            .collect()
            .map{ ids -> [ids] }  // wrap to prevent spreading during combine
            .combine(isolates_needing_update)
            .map{ all_ids, update_ids ->
                def all_list    = [all_ids].flatten()
                def update_list = [update_ids].flatten()
                def skipped  = all_list - update_list
                def will_run = all_list.intersect(update_list)

                println "${blue}=======================================================${reset}"
                println "${blue}  UPDATE_PHOENIX: The following isolate(s) WILL be${reset}"
                println "${blue}  processed through the update steps:${reset}"
                if (will_run) {
                    will_run.each { id -> println "${blue}    - ${id}${reset}" }
                } else {
                    println "${blue}    (none — all isolates are already current)${reset}"
                }
                println "${blue}=======================================================${reset}"

                println "${orange}=======================================================${reset}"
                println "${orange}  UPDATE_PHOENIX: The following isolate(s) already have${reset}"
                println "${orange}  the newest database files and will be skipped:${reset}"
                if (skipped) {
                    skipped.each { id -> println "${orange}    - ${id}${reset}" }
                } else {
                    println "${orange}    (none — all isolates require updating)${reset}"
                }
                println "${orange}=======================================================${reset}"
                return "printed"
            }
            .view{ "" }

        ///// RUN SHIGAPASS IF IT WASN'T BEFORE AND IS CORRECT TAXA /////

        // First, create a set of isolate IDs that already have shigapass files --> we will use this to filter and only run samples that don't have shigapass files and need them
        existing_shigapass_ids = CREATE_INPUT_CHANNELS.out.shigapass.map{ meta, shigapass_file -> meta.id}.ifEmpty(["none"]).toList() //[[id:], []]

        scaffolds_and_taxa_ch = CREATE_INPUT_CHANNELS.out.taxonomy.map{it -> get_taxa(it)}.filter{it, meta, taxonomy -> it.contains("Escherichia") || it.contains("Shigella")}.map{get_taxa_output, meta, taxonomy -> [[id:meta.id, project_id:meta.project_id], taxonomy ]}
                    .join(scaffolds_for_update_ch.map{                       meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}, by: [[0][0],[0][1]])
                    .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.map{ meta, fairy_outcome ->
                            // Read the content of the fairy_outcome file to check for FAILED
                            def content = file(fairy_outcome.toString()).text
                            def passed = !content.contains("FAILED")
                            [[id:meta.id, project_id:meta.project_id], [fairy_outcome, passed]]}, by: [[0][0],[0][1]])
                            .filter{ meta, taxonomy, filtered_scaffolds, fairy_data -> 
                                fairy_data[1] == true // Keep only entries where passed is true (no FAILED found)
                            }.map{   meta, taxonomy, filtered_scaffolds, fairy_data -> [meta, taxonomy, filtered_scaffolds] }.combine(existing_shigapass_ids.map{ [it] })
                            .filter{ meta, taxonomy, filtered_scaffolds, existing_shigapass_ids -> existing_shigapass_ids == 'none' || !existing_shigapass_ids.contains(meta.id)}
                            .map{    meta, taxonomy, filtered_scaffolds, existing_shigapass_ids -> [meta, taxonomy, filtered_scaffolds ] }// Add the filter to exclude isolates that already have shigapass files

        // Get ID from ShigaPass
        SHIGAPASS (
            scaffolds_and_taxa_ch, params.shigapass_database
        )
        ch_versions = ch_versions.mix(SHIGAPASS.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        checking_taxa_ch = CREATE_INPUT_CHANNELS.out.ani_best_hit.map{meta, ani_best_hit -> [[id:meta.id, project_id:meta.project_id], ani_best_hit]}
            .join(CREATE_INPUT_CHANNELS.out.ani.map{                  meta, ani          -> [[id:meta.id, project_id:meta.project_id], ani ]},      by: [[0][0],[0][1]])
            .join(SHIGAPASS.out.summary.map{                          meta, summary      -> [[id:meta.id, project_id:meta.project_id], summary ]},  by: [[0][0],[0][1]])
            .join(CREATE_INPUT_CHANNELS.out.taxonomy.map{             meta, taxonomy     -> [[id:meta.id, project_id:meta.project_id], taxonomy ]}, by: [[0][0],[0][1]])

        // check shigapass and correct fastani taxa if its wrong
        CHECK_SHIGAPASS_TAXA (
            checking_taxa_ch
        )
        ch_versions = ch_versions.mix(CHECK_SHIGAPASS_TAXA.out.versions)

        // If there was a change because of  assembly stats based on meta.id
        assembly_ratios_ch = CHECK_SHIGAPASS_TAXA.out.tax_file.map{meta, tax_file     -> [[id:meta.id, project_id:meta.project_id], tax_file]}\
            .join(CREATE_INPUT_CHANNELS.out.quast_report.map{      meta, quast_report -> [[id:meta.id, project_id:meta.project_id], quast_report]}, by: [[0][0],[0][1]])

        // Calculating the assembly ratio and gather GC% stats
        CALCULATE_ASSEMBLY_RATIO (
            assembly_ratios_ch, params.ncbi_assembly_stats
        )
        ch_versions = ch_versions.mix(CALCULATE_ASSEMBLY_RATIO.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file - written differently as fairy files might be variable in the number of lines in the file
        filtered_scaffolds_ch = scaffolds_for_update_ch.map{ meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}
                                .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.map{ meta, fairy_outcome ->
                                    // Read the content of the fairy_outcome file to check for FAILED
                                    def content = file(fairy_outcome.toString()).text
                                    def passed = !content.contains("FAILED")
                                    [[id:meta.id, project_id:meta.project_id], [fairy_outcome, passed]]}, by: [[0][0],[0][1]])
                                .filter { meta, filtered_scaffolds, fairy_data -> 
                                    fairy_data[1] // Keep only entries where passed is true (no FAILED found)
                                    }.map{ meta, filtered_scaffolds, fairy_data -> [meta, filtered_scaffolds] }

        // Running gamma to identify AR genes in scaffolds
        GAMMA_AR (
            filtered_scaffolds_ch, params.ardb
        )
        ch_versions = ch_versions.mix(GAMMA_AR.out.versions)

        // Running gamma to identify HV genes in scaffolds
        GAMMA_HV (
            filtered_scaffolds_ch, params.hvgamdb
        )
        ch_versions = ch_versions.mix(GAMMA_HV.out.versions)
        
        // Running gamma to identify plasmid replicons in scaffolds
        GAMMA_PF (
            filtered_scaffolds_ch, params.gamdbpf
        )
        ch_versions = ch_versions.mix(GAMMA_PF.out.versions)

        // combing fastp_trimd information with fairy check of reads to confirm there are reads after filtering
        trimd_reads_file_integrity_ch = reads_for_update_ch
            .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.map{ meta, fairy_outcome ->
                def content = file(fairy_outcome.toString()).text
                def passed = !content.contains("FAILED")
                [[id:meta.id, project_id:meta.project_id], [fairy_outcome, passed]]}, by: [[0][0],[0][1]])
            .filter { meta, reads, fairy_data -> fairy_data[1] == true }
            .map{ meta, reads, fairy_data -> [meta, reads] }
            .join(CREATE_INPUT_CHANNELS.out.mode_type, by: [0])
            .filter{ meta, reads, et -> et.base == 'CDC_PHOENIX' }
            .map { meta, reads, mode_type -> [meta, reads] }

        // Identifying AR genes in trimmed reads - only runs for entry types that have reads
        SRST2_AR (
            trimd_reads_file_integrity_ch, "gene", params.ardb
        )
        ch_versions = ch_versions.mix(SRST2_AR.out.versions)

        DO_MLST (
            scaffolds_for_update_ch,
            CREATE_INPUT_CHANNELS.out.fairy_outcome,
            reads_for_update_ch
                .join(CREATE_INPUT_CHANNELS.out.mode_type, by: [0])
                .filter{ meta, reads, et -> et.base == 'CDC_PHOENIX' }
                .map{ meta, reads, mode_type -> [meta, reads] },
            CHECK_SHIGAPASS_TAXA.out.tax_file
                .concat(CREATE_INPUT_CHANNELS.out.taxonomy)
                .unique{ meta, file -> [meta.id, meta.project_id] }
                .combine(isolates_needing_update)
                .filter{ meta, tax_file, update_ids -> update_ids.contains(meta.id) }
                .map{    meta, tax_file, update_ids -> [meta, tax_file] },
            ASSET_CHECK.out.mlst_db,
            true, "update"
        )
        ch_versions = ch_versions.mix(DO_MLST.out.versions)

        // Create file that has the organism name to pass to AMRFinder
        GET_TAXA_FOR_AMRFINDER (
            CHECK_SHIGAPASS_TAXA.out.tax_file
                .concat(CREATE_INPUT_CHANNELS.out.taxonomy)
                .unique{ meta, file -> [meta.id, meta.project_id] }
                .combine(isolates_needing_update)
                .filter{ meta, tax_file, update_ids -> update_ids.contains(meta.id) }
                .map{    meta, tax_file, update_ids -> [meta, tax_file] }
        )
        ch_versions = ch_versions.mix(GET_TAXA_FOR_AMRFINDER.out.versions)

        // Combining taxa and scaffolds to run amrfinder and get the point mutations.
        amr_channel = scaffolds_for_update_ch.map{          meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}
        .join(GET_TAXA_FOR_AMRFINDER.out.amrfinder_taxa.splitCsv(strip:true).map{meta, amrfinder_taxa     -> [[id:meta.id, project_id:meta.project_id], amrfinder_taxa ]}, by: [[0][0],[0][1]])
        .join(CREATE_INPUT_CHANNELS.out.prokka_faa.map{                          meta, prokka_faa         -> [[id:meta.id, project_id:meta.project_id], prokka_faa ]},     by: [[0][0],[0][1]])
        .join(CREATE_INPUT_CHANNELS.out.prokka_gff.map{                          meta, prokka_gff         -> [[id:meta.id, project_id:meta.project_id], prokka_gff ]},     by: [[0][0],[0][1]])

        // Run AMRFinder
        AMRFINDERPLUS_RUN (
            amr_channel, params.amrfinder_db
        )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

        files_to_update_ch = scaffolds_for_update_ch
            .map { meta, scaffolds -> [[id:meta.id, project_id:meta.project_id], meta] }
            .join(CREATE_INPUT_CHANNELS.out.pipeline_info_isolate.map{ meta, file -> [[id:meta.id, project_id:meta.project_id], file] }, by: [[0][0],[0][1]])
            .map { key, meta, pipeline_info -> [meta, pipeline_info] }
            .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{ meta, dir -> [[id:meta.id, project_id:meta.project_id], file(dir)] }, by: [[0][0],[0][1]])
            .map{ meta, a, b -> 
                def pipeline_info = [a, b].find { it.toString().endsWith('.yml') }
                def dir = [a, b].find { !it.toString().endsWith('.yml') }
                [meta, dir, pipeline_info]
            }
            .join(CREATE_INPUT_CHANNELS.out.readme.map{                meta, readme -> [[id:meta.id, project_id:meta.project_id], readme]}, by: [[0][0],[0][1]], remainder: true)
            .filter { it -> it.size() == 4 }
            .map{ meta, dir, pipeline_info, readme -> [meta, dir, pipeline_info, readme ?: []] }
            .join(CREATE_INPUT_CHANNELS.out.gamma_ar.map{              meta, gamma_ar    -> [[id:meta.id, project_id:meta.project_id], gamma_ar]},    by: [[0][0],[0][1]])
            .join(GAMMA_AR.out.gamma.map{                              meta, gamma       -> [[id:meta.id, project_id:meta.project_id], gamma]},       by: [[0][0],[0][1]])
            .join(CREATE_INPUT_CHANNELS.out.ncbi_report.map{           meta, ncbi_report -> [[id:meta.id, project_id:meta.project_id], ncbi_report]}, by: [[0][0],[0][1]])
            .join(AMRFINDERPLUS_RUN.out.report.map{                    meta, report      -> [[id:meta.id, project_id:meta.project_id], report]},      by: [[0][0],[0][1]])
            .join(CREATE_INPUT_CHANNELS.out.taxonomy.map{              meta, taxonomy    -> [[id:meta.id, project_id:meta.project_id], taxonomy]},    by: [[0][0],[0][1]])
            .join(CHECK_SHIGAPASS_TAXA.out.edited_tax_file.map{        meta, tax_file    -> [[id:meta.id, project_id:meta.project_id], tax_file]},    by: [[0][0],[0][1]], remainder: true)
            .join(CREATE_INPUT_CHANNELS.out.gamma_pf.map{              meta, gamma_pf    -> [[id:meta.id, project_id:meta.project_id], gamma_pf]},    by: [[0][0],[0][1]])
            .join(GAMMA_PF.out.gamma.map{                              meta, g_pf        -> [[id:meta.id, project_id:meta.project_id], g_pf]},        by: [[0][0],[0][1]])
            .join(CREATE_INPUT_CHANNELS.out.gamma_hv.map{              meta, gamma_hv    -> [[id:meta.id, project_id:meta.project_id], gamma_hv]},    by: [[0][0],[0][1]])
            .join(GAMMA_HV.out.gamma.map{                              meta, g_hv        -> [[id:meta.id, project_id:meta.project_id], g_hv]},        by: [[0][0],[0][1]])
            .join(CREATE_INPUT_CHANNELS.out.update_pipeline_info_isolate
                .combine(isolates_needing_update)
                .filter { meta, software, update_ids -> 
                    def flat = [update_ids].flatten().findAll { it != null && it != "" }
                    flat.contains(meta.id)
                }
                .map { meta, software, update_ids -> [[id:meta.id, project_id:meta.project_id], software] }
            , by: [[0][0],[0][1]], remainder: true)
            .map { meta, dir, pipeline_info, readme, g_ar, g, ncbi, rpt, tax, t_file, g_pf, gp, g_hv, gh, software ->
                return [
                    meta, dir, pipeline_info, readme, g_ar, g, ncbi, rpt, tax, t_file ?: [], g_pf, gp, g_hv, gh, software ?: []
                ]
            }

        files_to_update_local = files_to_update_ch

        CREATE_AND_UPDATE_README (
            files_to_update_local,
            wf_version,
            params.custom_mlstdb,
            params.ardb,
            params.gamdbpf,
            params.hvgamdb,
            params.amrfinder_db
        )
        ch_versions = ch_versions.mix(CREATE_AND_UPDATE_README.out.versions)

        // Create scaffold placeholder channel - one empty entry per scaffold sample
        scaffold_meta_ch = CREATE_INPUT_CHANNELS.out.mode_type
            .combine(isolates_needing_update)
            .filter{ meta, et, update_ids -> update_ids.contains(meta.id) }
            .filter{ meta, et, update_ids -> et.base == 'SCAFFOLDS' || et.base == 'CDC_SCAFFOLDS' }
            .map{    meta, et, update_ids -> [meta, []] }

        // Filter mode_type to only samples that need updating
        ch_mode_type = filter_to_update_ids(
            CREATE_INPUT_CHANNELS.out.mode_type,
            isolates_needing_update
        )

        // Filter concat channels to only samples that need updating
        ch_taxonomy = filter_to_update_ids(
            merge_unique(
                CHECK_SHIGAPASS_TAXA.out.tax_file,
                CREATE_INPUT_CHANNELS.out.taxonomy
            ),
            isolates_needing_update
        )

        ch_ani_best_hit = filter_to_update_ids(
            merge_unique(
                CHECK_SHIGAPASS_TAXA.out.ani_best_hit,
                CREATE_INPUT_CHANNELS.out.ani_best_hit
            ),
            isolates_needing_update
        )

        GENERATE_PIPELINE_STATS_WF (
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.raw_stats, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.fastp_total_qc, isolates_needing_update),
            filter_to_update_ids(SRST2_AR.out.fullgene_results.concat(CREATE_INPUT_CHANNELS.out.srst2_ar).unique{ meta, file -> [meta.id, meta.project_id] }, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_trimd_report, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_trimd_krona, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_trimd_bh_summary, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.renamed_scaffolds, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.filtered_scaffolds, isolates_needing_update),
            filter_to_update_ids(DO_MLST.out.checked_MLSTs, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.gamma_hv, isolates_needing_update),
            filter_to_update_ids(GAMMA_AR.out.gamma, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.gamma_pf, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.quast_report, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.busco_short_summary, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_asmbld_report, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_asmbld_krona, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_asmbld_bh_summary, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_wtasmbld_report, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_wtasmbld_krona, isolates_needing_update),
            filter_to_update_ids(CREATE_INPUT_CHANNELS.out.k2_wtasmbld_bh_summary, isolates_needing_update),
            filter_to_update_ids(CHECK_SHIGAPASS_TAXA.out.tax_file.concat(CREATE_INPUT_CHANNELS.out.taxonomy).unique{ meta, file -> [meta.id, meta.project_id] }, isolates_needing_update),
            filter_to_update_ids(CHECK_SHIGAPASS_TAXA.out.ani_best_hit.concat(CREATE_INPUT_CHANNELS.out.ani_best_hit).unique{ meta, file -> [meta.id, meta.project_id] }, isolates_needing_update),
            filter_to_update_ids(CALCULATE_ASSEMBLY_RATIO.out.ratio.concat(CREATE_INPUT_CHANNELS.out.assembly_ratio).unique{ meta, file -> [meta.id, meta.project_id] }, isolates_needing_update),
            filter_to_update_ids(AMRFINDERPLUS_RUN.out.mutation_report, isolates_needing_update),
            filter_to_update_ids(CALCULATE_ASSEMBLY_RATIO.out.gc_content.concat(CREATE_INPUT_CHANNELS.out.gc_content).unique{ meta, file -> [meta.id, meta.project_id] }, isolates_needing_update),
            ch_mode_type
        )
        ch_versions = ch_versions.mix(GENERATE_PIPELINE_STATS_WF.out.versions)

        //get back project_id info
        pipeline_stats_ch = GENERATE_PIPELINE_STATS_WF.out.pipeline_stats.concat(CREATE_INPUT_CHANNELS.out.synopsis).unique{ meta, file -> [meta.id, meta.project_id] }
                                    .join(GAMMA_AR.out.gamma.map{ meta, gamma -> [[id:meta.id], gamma, meta.project_id]}, by: [0]).map{ meta, pipeline_stats, gamma, project_id -> [[id:meta.id, project_id:project_id], pipeline_stats]}


        // Combining output based on meta.id to create summary by sample
        // Anchor: scaffolds_for_update_ch — present for ALL sample types (scaffolds + reads)
        scaffolds_local = scaffolds_for_update_ch
            .map{ meta, scaffolds -> [[id:meta.id, project_id:meta.project_id], scaffolds] }
        line_step2 = scaffolds_local.join(DO_MLST.out.checked_MLSTs
                .map{ meta, checked_MLSTs -> [[id:meta.id, project_id:meta.project_id], checked_MLSTs] }, by: [[0][0],[0][1]])
        line_step3 = line_step2.join(CREATE_INPUT_CHANNELS.out.gamma_hv
                .map{ meta, gamma_hv -> [[id:meta.id, project_id:meta.project_id], gamma_hv] }, by: [[0][0],[0][1]])
        line_step4 = line_step3.join(GAMMA_AR.out.gamma
                .map{ meta, gamma -> [[id:meta.id, project_id:meta.project_id], gamma] }, by: [[0][0],[0][1]])
        line_step5 = line_step4.join(CREATE_INPUT_CHANNELS.out.gamma_pf
                .map{ meta, gamma_pf -> [[id:meta.id, project_id:meta.project_id], gamma_pf] }, by: [[0][0],[0][1]])
        line_step6 = line_step5.join(CREATE_INPUT_CHANNELS.out.quast_report
                .map{ meta, quast_report -> [[id:meta.id, project_id:meta.project_id], quast_report] }, by: [[0][0],[0][1]])
        line_step7 = line_step6.join(CALCULATE_ASSEMBLY_RATIO.out.ratio
                .concat(CREATE_INPUT_CHANNELS.out.assembly_ratio)
                .unique{ meta, file -> [meta.id, meta.project_id] }
                .map{ meta, assembly_ratio -> [[id:meta.id, project_id:meta.project_id], assembly_ratio] }, by: [[0][0],[0][1]])
        line_step8 = line_step7.join(GENERATE_PIPELINE_STATS_WF.out.pipeline_stats
                .concat(CREATE_INPUT_CHANNELS.out.synopsis)
                .unique{ meta, file -> [meta.id, meta.project_id] }
                .map{ meta, synopsis -> [[id:meta.id, project_id:meta.project_id], synopsis] }, by: [[0][0],[0][1]])
        line_step9 = line_step8.join(CHECK_SHIGAPASS_TAXA.out.tax_file
                .concat(CREATE_INPUT_CHANNELS.out.taxonomy)
                .unique{ meta, file -> [meta.id, meta.project_id] }
                .map{ meta, taxonomy -> [[id:meta.id, project_id:meta.project_id], taxonomy] }, by: [[0][0],[0][1]])
        line_step10 = line_step9.join(CREATE_INPUT_CHANNELS.out.k2_wtasmbld_bh_summary
                .map{ meta, k2_wtasmbld_bh_summary -> [[id:meta.id, project_id:meta.project_id], k2_wtasmbld_bh_summary] }, by: [[0][0],[0][1]])
        line_step11 = line_step10.join(AMRFINDERPLUS_RUN.out.report
                .map{ meta, report -> [[id:meta.id, project_id:meta.project_id], report] }, by: [[0][0],[0][1]])
        line_step12 = line_step11.join(CHECK_SHIGAPASS_TAXA.out.ani_best_hit
                .concat(CREATE_INPUT_CHANNELS.out.ani_best_hit)
                .unique{ meta, file -> [meta.id, meta.project_id] }
                .map{ meta, ani_best_hit -> [[id:meta.id, project_id:meta.project_id], ani_best_hit] }, by: [[0][0],[0][1]])
        line_step13 = line_step12.join(CREATE_INPUT_CHANNELS.out.pipeline_info_isolate
                .filter{ meta, pipeline_info -> pipeline_info != null }
                .map{ meta, pipeline_info ->
                    def content = pipeline_info.text
                    def version = "unknown"
                    if      (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/) { version = (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/)[0][1].trim() }
                    else if (content =~ /(?m)version:\s*(.+)/)          { version = (content =~ /(?m)version:\s*(.+)/)[0][1].trim() }
                    else if (content =~ /(?m)phoenix:\s*(.+)/)          { version = (content =~ /(?m)phoenix:\s*(.+)/)[0][1].trim() }
                    [[id:meta.id, project_id:meta.project_id], version] }, by: [[0][0],[0][1]])

        // ── Null-safe late joins for read level outputs ──────────────────────────
        line_step14 = line_step13.join(CREATE_INPUT_CHANNELS.out.fastp_total_qc
                .map{ meta, f -> [[id:meta.id, project_id:meta.project_id], f] }, by: [[0][0],[0][1]], remainder: true)
            .filter { it -> it[1] != null }  // drop right-side remainders
            .map{ it -> it[-1] == null ? it[0..-2] + [[]] : it }

        line_step15 = line_step14.join(CREATE_INPUT_CHANNELS.out.k2_trimd_bh_summary
                .map{ meta, k2_trimd_bh_summary -> [[id:meta.id, project_id:meta.project_id], k2_trimd_bh_summary] }, by: [[0][0],[0][1]], remainder: true)
            .filter { it -> it[1] != null }  // drop right-side remainders
            .map{ it -> it[-1] == null ? it[0..-2] + [[]] : it }


       line_summary_ch = line_step15
            .map{ meta, scaffolds, checked_MLSTs, gamma_hv, gamma, gamma_pf, quast_report,
                        assembly_ratio, synopsis, taxonomy, k2_wtasmbld_bh_summary, report,
                        ani_best_hit, version, fastp_total_qc, k2_trimd_bh_summary ->
                [meta, fastp_total_qc, checked_MLSTs, gamma_hv, gamma, gamma_pf, quast_report,
                 assembly_ratio, synopsis, taxonomy, k2_trimd_bh_summary, k2_wtasmbld_bh_summary,
                 report, ani_best_hit, version]
            }

/*        line_summary_ch.view { it ->
            def meta = it[0]
            return """
            ====================================================
            SURVIVOR FOUND: ${meta.id}
            Total elements in tuple: ${it.size()}
            Last element check: ${it[-1]}
            ====================================================
            """.stripIndent()
        }
        line_summary_ch.ifEmpty { "!!! LOG ALERT: No samples survived the line_step joins !!!" }.view()
*/
        // First, check if SHIGAPASS.out.summary is empty and create appropriate channel
        shigapass_ch = SHIGAPASS.out.summary.mix(CREATE_INPUT_CHANNELS.out.shigapass)

        // Extract [id, project_id] pairs for all isolates in line_summary_ch channel for joining // If shigapass_files is null (no match), replace with empty list
        all_id_pairs_ch = line_summary_ch.join(shigapass_ch, by: [[0][0],[0][1]], remainder: true).filter{ it[1] != null }
            .map{ meta, fastp_total_qc, checked_MLSTs, gamma_hv, gamma, gamma_pf, quast_report, assembly_ratio, synopsis, taxonomy, k2_trimd_bh_summary, k2_wtasmbld_bh_summary, report, ani_best_hit, old_version, shigapass_file ->
            def sp = (shigapass_file == null || shigapass_file == []) ? [] : shigapass_file
            return [meta, fastp_total_qc, checked_MLSTs, gamma_hv, gamma, gamma_pf, quast_report, assembly_ratio, synopsis, taxonomy, k2_trimd_bh_summary, k2_wtasmbld_bh_summary, report, ani_best_hit, sp, old_version]
            }
        // rename to line_summary_ch --> just to keep the naming consistent with the original workflow
        line_summary_ch_final = all_id_pairs_ch

        // Generate summary per sample
        CREATE_SUMMARY_LINE (
            line_summary_ch_final, false, wf_version
        )
        ch_versions = ch_versions.mix(CREATE_SUMMARY_LINE.out.versions)

        // Extract the list of project folders from input_dir channel, combine with the actual path that is needed by GATHER_SUMMARY_LINES
        project_ids = CREATE_INPUT_CHANNELS.out.directory_ch.map{ it[1] }.unique().collect()
            .map { dirs -> dirs.collectEntries { dir -> def dirName = dir.tokenize('/').last()
            return [dirName, dir]}}

        // Implement fallback mechanism: if CREATE_SUMMARY_LINE produces null summaryline, use CREATE_INPUT_CHANNELS.out.line_summary
        fallback_line_summary_ch = CREATE_INPUT_CHANNELS.out.line_summary.map{meta, line_summary -> [[id:meta.id, project_id:meta.project_id], line_summary]}
        
        // First, join CREATE_SUMMARY_LINE with gene_results and identify entries with null summaryline
        primary_summaries_ch = CREATE_SUMMARY_LINE.out.line_summary
            .join(SRST2_AR.out.fullgene_results.map{meta, fullgene_results -> [[id:meta.id, project_id:meta.project_id], fullgene_results]}, by: [[0][0],[0][1]], remainder: true) // add in SRST2 results if they exist to make the pipeline wait for it to finish
            .map{meta, summaryline, fullgene_results -> [meta, summaryline, fullgene_results]}

        // Bring in existing summarylines for samples that were skipped (already current)
        existing_summaries_for_skipped_ch = CREATE_INPUT_CHANNELS.out.line_summary
            .combine(isolates_needing_update)
            .filter{ meta, summary, update_ids -> !update_ids.contains(meta.id) }
            .map{    meta, summary, update_ids -> [meta.project_id, summary] }
        
        // Split into successful and failed summaries based on whether summaryline is null
        successful_summaries_ch = primary_summaries_ch.filter{ meta, summaryline, fullgene_results -> summaryline != null }.map{meta, summaryline, fullgene_results -> [meta.project_id, summaryline]}
        failed_summaries_ch = primary_summaries_ch.filter{ meta, summaryline, fullgene_results -> summaryline == null }.map{meta, summaryline, fullgene_results -> [[id:meta.id, project_id:meta.project_id], fullgene_results]}

        // Implement fallback mechanism: if CREATE_SUMMARY_LINE produces null summaryline, use CREATE_INPUT_CHANNELS.out.line_summary
        // Join failed summaries with fallback line_summary
        fallback_summaries_ch = failed_summaries_ch.join(CREATE_INPUT_CHANNELS.out.line_summary.map{meta, line_summary -> [[id:meta.id, project_id:meta.project_id], line_summary]}, by: [[0][0],[0][1]])
            .map{meta, fullgene_results, line_summary -> [meta.project_id, line_summary]}

        // Combine successful and fallback summaries
        summaries_ch = successful_summaries_ch.mix(fallback_summaries_ch)
            .join(CREATE_INPUT_CHANNELS.out.mode_type
            .map{ meta, et -> [meta.project_id, et] }, by: [0])
            .groupTuple(by: 0)
            .map { group -> 
                def meta = [:]
                def (id, summary_lines, mode_types) = group
                id2 = id.split('/')[-1]
                meta.project_id = id.split('/')[-1]
                meta.full_project_id = id
                def full_project_id = project_ids.value[id2]
                def first_type = mode_types.flatten().first().base
                def busco_boolean = first_type == 'CDC_PHOENIX' || first_type == 'CDC_SCAFFOLDS'
                return [meta, summary_lines, full_project_id, busco_boolean]
            }

        //add in the pipeline info to get the version
        gathered_summaries_ch = summaries_ch.join(CREATE_INPUT_CHANNELS.out.pipeline_info.map{meta, pipeline_info -> [[project_id:meta.project_id.toString().split('/')[-1], full_project_id:meta.project_id], pipeline_info]}, by: [0])

        // For cases where the user is running isolates from multiple directories and wants the summary files output to their original project_id folders. 
        multiple_directories_ch = CREATE_INPUT_CHANNELS.out.directory_ch.map{meta, dir -> [dir]}.collect().map{ files -> files.unique().size() > 1 }
        // Conditionally group by project based on if there are multiple directories
        branched_collected_summaries_ch = gathered_summaries_ch
            .combine(multiple_directories_ch)
            .branch { meta, summary_line, dir, busco_boolean, pipeline_info, is_multiple ->
                multi_dir: is_multiple == true
                    return [meta, summary_line, dir, busco_boolean, pipeline_info]
                single_dir: is_multiple == false
                    return [meta, summary_line, dir, busco_boolean, pipeline_info] }

        def collected_summaries_ch
        if (params.outdir != "${launchDir}/phx_output") {
            //if there are multiple dirs and an --outdir given then the will need to all be combined into one Phoenix_Summary.tsv file in the outdir 
            collected_summaries_multi_ch = branched_collected_summaries_ch.multi_dir.collect().map{ items ->
                                def grouped_items = items.collate(5)
                                def summary_lines = []
                                def dirs = []
                                def busco_booleans = []
                                def pipeline_infos = []
                            grouped_items.each{ item ->
                                summary_lines.addAll(item[1])  // Use addAll to flatten the list
                                dirs.add(item[2])
                                busco_booleans.add(item[3])
                                pipeline_infos.add(item[4])}
                            return [[], summary_lines, dirs, busco_booleans.any{ it == true }, pipeline_infos]}

            //for single dirs
            collected_summaries_ch = branched_collected_summaries_ch.single_dir.mix(collected_summaries_multi_ch)  // take first dir/busco since they should be the same within a project

        } else { // **************************** this needs to be checked for multiple dirs  ****************************

            // Group by project and collect files within each project -- multi dir
            collected_summaries_multi_ch = branched_collected_summaries_ch.multi_dir
                .groupTuple(by: 0)  // group by meta (project_id)
                .map{ meta, summary_lines, dirs, busco_booleans, pipeline_info -> 
                    //def selected_summary_lines = summary_lines.size() > 1 ? summary_lines[0] : summary_lines
                    [meta, summary_lines, dirs.first().toString(), busco_booleans.any{ it == true }, pipeline_info]}  // take first dir/busco since they should be the same within a project

            // Group by project and collect files within each project -- single dir
            collected_summaries_ch = branched_collected_summaries_ch.single_dir.mix(collected_summaries_multi_ch)
        }

        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            collected_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> meta},
            collected_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> summary_lines.flatten()},
            collected_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> full_project_id},
            collected_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> busco_boolean},
            collected_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> pipeline_info}.map { file -> (file.text =~ /(?m)cdcgov\/phoenix: (.+)/)[0][1].trim() } // Extract the version from the pipeline_info file
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        //this means that files need to be directed to the --outdir so we need to update the dir in the channel 
        outdir_full_path = Channel.fromPath(params.outdir, type: 'dir') // get the full path to the outdir, by not using "relative: true"
        summaries_with_outdir_ch = summaries_ch.combine(outdir_full_path.toList()).map{meta, summary_lines, full_project_id, busco_boolean, outdir -> [meta, summary_lines, outdir, busco_boolean] }

        // Now we need to check if --centar was passed when the samples were run previously. // "ifEmpty()" branch executes if no files match 
        centar_var = CREATE_INPUT_CHANNELS.out.mode_type
            .map{ meta, et -> et.centar }
            .filter{ it == true }
            .ifEmpty( Channel.of(false) )
            .first()

        srst2_for_griphin = SRST2_AR.out.fullgene_results
            .join(reads_for_update_ch.map { meta, reads -> [meta.id, true] }, by: [0])
            .map { meta, f, flag -> [meta, f] }

        version_per_project_ch = collected_summaries_ch
            .map { meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> 
                def content = pipeline_info.text
                def version = "unknown"
                
                // Try multiple patterns
                if (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/) {
                    version = (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/)[0][1].trim()
                } else if (content =~ /(?m)version:\s*(.+)/) {
                    version = (content =~ /(?m)version:\s*(.+)/)[0][1].trim()
                } else if (content =~ /(?m)phoenix:\s*(.+)/) {
                    version = (content =~ /phoenix:\s*(.+)/)[0][1].trim()
                }
                
               // Try matching the key structure - use full_project_id if that's what griphin_inputs_ch uses
                [[project_id: full_project_id], version]  // Changed from meta.project_id to full_project_id
            }

        //create GRiPHin report channel
        griphin_inputs_ch = empty_ch
            .mix(
                CREATE_INPUT_CHANNELS.out.fastp_total_qc,
                CREATE_INPUT_CHANNELS.out.raw_stats,
                CREATE_INPUT_CHANNELS.out.k2_trimd_bh_summary,
                CREATE_INPUT_CHANNELS.out.k2_trimd_report,
                CREATE_INPUT_CHANNELS.out.k2_wtasmbld_bh_summary,
                CREATE_INPUT_CHANNELS.out.k2_wtasmbld_report,
                CREATE_INPUT_CHANNELS.out.quast_report,
                CREATE_INPUT_CHANNELS.out.fairy_outcome,
                DO_MLST.out.checked_MLSTs,
                CHECK_SHIGAPASS_TAXA.out.tax_file.concat(CREATE_INPUT_CHANNELS.out.taxonomy).unique{ meta, file -> [meta.id, meta.project_id] },
                CALCULATE_ASSEMBLY_RATIO.out.ratio.concat(CREATE_INPUT_CHANNELS.out.assembly_ratio).unique{ meta, file -> [meta.id, meta.project_id] },
                CALCULATE_ASSEMBLY_RATIO.out.gc_content.concat(CREATE_INPUT_CHANNELS.out.gc_content).unique{ meta, file -> [meta.id, meta.project_id] },
                GAMMA_AR.out.gamma,
                GAMMA_PF.out.gamma,
                GAMMA_HV.out.gamma,
                CHECK_SHIGAPASS_TAXA.out.ani_best_hit.concat(CREATE_INPUT_CHANNELS.out.ani_best_hit).unique{ meta, file -> [meta.id, meta.project_id] },
                pipeline_stats_ch.concat(CREATE_INPUT_CHANNELS.out.synopsis).unique{ meta, file -> [meta.id, meta.project_id] },
                SHIGAPASS.out.summary.concat(CREATE_INPUT_CHANNELS.out.shigapass).unique{ meta, file -> [meta.id, meta.project_id] },
                CREATE_INPUT_CHANNELS.out.busco_short_summary,
                srst2_for_griphin
            )

        def software_versions_ch =empty_ch
        if (params.indir != null) {
            // Use the existing collected_summaries_ch - no need to recalculate
            shigapass_var = CHECK_SHIGAPASS_TAXA.out.tax_file.concat(CREATE_INPUT_CHANNELS.out.taxonomy)
                .unique{ meta, file -> [meta.id, meta.project_id] }
                .map{it -> get_only_taxa(it)}.collect().flatten()
                .count{ it -> it.contains("Escherichia") || it.contains("Shigella")}
                .collect().sum().map{ it -> it[0] > 0 }

            busco_boolean = collected_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> busco_boolean}
            
            def outdir_full_path2
            if (params.outdir != "${launchDir}/phx_output") {
                outdir_full_path2 = Channel.fromPath(params.outdir, type: 'dir')
            } else {
                outdir_full_path2 = Channel.fromPath(params.indir, type: 'dir')
            }

            old_versions_collected = collected_summaries_ch
                .map { meta, files, path, busco, info -> 
                    def version = "unknown"
                    if (info instanceof Path && info.exists()) {
                        def content = info.text
                        if (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/) {
                            version = (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/)[0][1].trim()
                        }
                    }
                    return version
                }
 
            processed_griphin_inputs = griphin_inputs_ch.filter { meta, file -> file != null }
                .groupTuple()
                .map { meta, files ->
                    def is_scaffold = !files.any { it.name.contains("fastp") }
                    [
                        meta: [ 
                            id: "${meta.id}", 
                            filenames: files.collect { it.getName() },
                            is_scaffold: is_scaffold
                        ],
                        files: files 
                    ] 
                }

            // Sets inferred mode for GRiPHin based on whether the input is from scaffolds or reads, determined by checking if any of the files for that sample contain "fastp" in the name (indicating read-level data). If no files contain "fastp", it is assumed to be scaffold-level data and inferred mode is set to "SCAFFOLD_INFERRED". If any file contains "fastp", inferred mode is set to "READS". This allows GRiPHin to adjust its processing based on the type of input data without requiring an explicit parameter from the user.
            indir_inferred_mode = CREATE_INPUT_CHANNELS.out.mode_type
                .map{ meta, et -> et.base }
                .collect()
                .map{ bases ->
                    def has_reads = bases.any{ it == 'CDC_PHOENIX' || it == 'PHOENIX' }
                    has_reads ? "READS" : "SCAFFOLD_INFERRED"
                }

            GRIPHIN_PUBLISH (
                params.ardb,
                CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                processed_griphin_inputs.map { it.meta }.collect(),
                processed_griphin_inputs.map { it.files }.collect(),
                outdir_full_path2,
                workflow.manifest.version,
                params.coverage,
                busco_boolean,
                shigapass_var, 
                false, 
                params.bldb, 
                false, 
                false, 
                [],
                old_versions_collected,
                indir_inferred_mode
            )
            
            ch_versions = ch_versions.mix(GRIPHIN_PUBLISH.out.versions)
            griphin_tsv_report = GRIPHIN_PUBLISH.out.griphin_tsv_report
            griphin_report = GRIPHIN_PUBLISH.out.griphin_report

            ch_collated_versions = ch_versions
                .unique()
                .collectFile(name: 'collated_versions.yml')
                .ifEmpty { 
                    log.info ">>> ch_collated_versions was empty, using fallback"
                    Channel.fromPath("${workflow.projectDir}/assets/nf-core_version.yml")
                }
                
                software_versions_ch = CREATE_INPUT_CHANNELS.out.directory_ch
                    .combine(isolates_needing_update)
                    .filter { meta, dir, update_ids -> update_ids.contains(meta.id) }
                .map{ meta, directory_ch, update_ids -> 
                    [[project_id: meta.project_id.toString().split('/')[-1].replace("]", ""), 
                    full_project_id: directory_ch], 
                    directory_ch]
                }
                .unique()
                .combine(ch_collated_versions)
                .map { meta, dir, v_file -> [meta, v_file] }

            META_CUSTOM_DUMPSOFTWAREVERSIONS ( software_versions_ch )

        } else { // for --input
            // 1. BRIDGE: Define shigapass check (MOVED HERE TO BE ACCESSIBLE TO BOTH BRANCHES)
            ch_shigapass_check = CHECK_SHIGAPASS_TAXA.out.tax_file
                .mix(CREATE_INPUT_CHANNELS.out.taxonomy)
                .map{ get_only_taxa(it) }
                .collect()
                .map{ it.flatten().any { tax -> tax.contains("Escherichia") || tax.contains("Shigella") } }
                .ifEmpty(false)

            if (params.outdir == "${launchDir}/phx_output") {
                // BRANCH 3: --input WITHOUT --outdir (GRIPHIN_NO_PUBLISH per-project)
               
                // Extract old_version PER PROJECT from collected_summaries_ch
                ch_old_versions_per_project = collected_summaries_ch
                    .map { meta, files, path, busco, info -> 
                        def pid = meta.project_id.toString().split('/')[-1].replace("]", "").trim()

                        def infoFile = info instanceof Collection ? info[0] : info

                        def version = "unknown"
                        if (infoFile instanceof Path && infoFile.exists()) {
                            def content = infoFile.text
                            if (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/) {
                                version = (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/)[0][1].trim()
                            }
                        }
                        return [pid, version]
                    }

                // Determine inferred_mode PER PROJECT
                ch_inferred_mode_per_project = CREATE_INPUT_CHANNELS.out.mode_type
                    .map { meta, mode_type ->
                        def pid = meta.project_id.toString().split('/')[-1].replace("]", "").trim()
                        def mode = (mode_type.base in ['CDC_PHOENIX', 'PHOENIX']) ? 'READS' : 'SCAFFOLD_INFERRED'
                        [pid, mode]
                    }
                    .groupTuple(by: 0)
                    .map { pid, modes ->
                        [pid, modes.flatten().contains('READS') ? 'READS' : 'SCAFFOLD_INFERRED']
                    }

                ch_per_project = griphin_inputs_ch
                    .filter{ meta, file -> file != null }  // same filter here
                    .map { meta, file ->
                        def pid = meta.project_id.toString().split('/')[-1].replace("]", "").trim()
                        return [ pid, meta, file ]
                    }
                    .groupTuple(by: 0)
                    .map { pid, metas, files ->
                        def clean_files = files.flatten().findAll{ it != null }
                        def sample_ids = metas.unique{ it.id }*.id
                        def meta_list = sample_ids.collect { sid ->
                            [ id: sid, filenames: clean_files.findAll { it.name.contains(sid) }.collect { it.name } ]
                        }
                        return [ pid, meta_list, clean_files ]
                    }

                // Get busco per project
                ch_busco_per_project = collected_summaries_ch
                    .map { meta, files, path, busco, info -> 
                        def pid = meta.project_id.toString().split('/')[-1].replace("]", "").trim()
                        return [ pid, busco, path ] 
                    }
/*
                // 1) Show exact keys for every project-level channel
                ch_per_project
                    .map { pid, meta_list, clean_files -> "|${pid}|" }
                    .view { "PER_PROJECT KEY: ${it}" }

                ch_busco_per_project
                    .map { pid, busco, path -> "|${pid}|" }
                    .view { "BUSCO KEY: ${it}" }

                ch_old_versions_per_project
                    .map { pid, old_version -> "|${pid}|" }
                    .view { "OLD_VERSION KEY: ${it}" }

                ch_inferred_mode_per_project
                    .map { pid, mode -> "|${pid}|" }
                    .view { "MODE KEY: ${it}" }

                // 2) Check whether any project-level channel has duplicate rows per project
                ch_busco_per_project
                    .groupTuple(by: 0)
                    .view { pid, buscos, paths ->
                        "BUSCO ROWS FOR ${pid}: count=${buscos.size()} buscos=${buscos}"
                    }

                ch_old_versions_per_project
                    .groupTuple(by: 0)
                    .view { pid, versions ->
                        "OLD_VERSION ROWS FOR ${pid}: count=${versions.size()} versions=${versions}"
                    }

                ch_per_project
                    .view { pid, meta_list, clean_files ->
                        "PER_PROJECT FULL ${pid}: meta_count=${meta_list.size()} file_count=${clean_files.size()} files=${clean_files*.name}"
                    }
*/
                // 3) Split joins into separate steps so you can see exactly where it stops
                j1 = ch_per_project
                    .join(ch_busco_per_project, by: 0)
                
                j2 = j1
                    .join(ch_old_versions_per_project, by: 0)
                
                j3 = j2
                    .join(ch_inferred_mode_per_project, by: 0)
                
                ch_combined = j3
                    .map { row ->
                        assert row.size() == 7 : "Unexpected ch_combined shape: ${row}"
                        def (pid, meta_list, files, busco, path, old_version, inferred_mode) = row
                        [meta_list, files, path, busco, old_version, inferred_mode]
                    }
                
                // 4) Name the GRIPHIN inputs explicitly and print them before the call
                ch_valid_samplesheet = CREATE_INPUT_CHANNELS.out.valid_samplesheet
                    .collect()
                    .map { it ?: [] }
 //
                ch_griphin_metas = ch_combined
                    .map { it[0] }
//                    .view { "GRIPHIN METAS: ${it}" }

                ch_griphin_files = ch_combined
                    .map { it[1] }
//                    .view { "GRIPHIN FILES: ${it*.name}" }

                ch_griphin_path = ch_combined
                    .map { it[2] }
//                    .view { "GRIPHIN PATH: ${it}" }

                ch_griphin_busco = ch_combined
                    .map { it[3] }
//                    .view { "GRIPHIN BUSCO: ${it}" }

                ch_griphin_old_version = ch_combined
                    .map { it[4] }
//                    .view { "GRIPHIN OLD_VERSION: ${it}" }

                ch_griphin_mode = ch_combined
                    .map { it[5] }
//                    .view { "GRIPHIN MODE: ${it}" }

                GRIPHIN_NO_PUBLISH (
                    params.ardb,
                    CREATE_INPUT_CHANNELS.out.valid_samplesheet.collect().ifEmpty([]),
                    ch_combined.map { it[0] },  // metas
                    ch_combined.map { it[1] },  // files
                    ch_combined.map { it[2] },  // path
                    workflow.manifest.version,
                    params.coverage,
                    ch_combined.map { it[3] },  // busco
                    ch_shigapass_check, 
                    false, 
                    params.bldb,
                    false, 
                    true,
                    [],
                    ch_combined.map { it[4] },  // old_version PER PROJECT
                    ch_combined.map { it[5] }   // inferred_mode PER PROJECT
                )
                ch_versions = ch_versions.mix(GRIPHIN_NO_PUBLISH.out.versions)

                // 3. Process the new reports and add metadata for joining
                griphin_reports_ch = GRIPHIN_NO_PUBLISH.out.griphin_report.collect().ifEmpty([[],[]]).flatten().collate(2)
                                    .map{ path_txt, griphin_report -> add_meta(path_txt, griphin_report) }

                // 4. Join old Excel, new report, and project directory
                griphins_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch
                    .map{ meta, excel -> [[project_id: meta.project_id.toString().split('/')[-1].replace("]", "")], excel] }.unique()
                    .join(griphin_reports_ch.map{ meta, report -> [[project_id: meta.project_id.toString().split('/')[-1].replace("]", "")], report] }, by: [0])
                    .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{ meta, dir -> [[project_id: meta.project_id.toString().split('/')[-1].replace("]", "")], dir] }.unique(), by: [0])
                    .map { list ->
                        return list
                    }
                    .map{ meta, old_excel, new_excel, directory -> [[project_id: meta.project_id, full_project_id: directory], old_excel, new_excel] }

                // 5. Merge old and new data
                UPDATE_GRIPHIN (
                    griphins_ch.map{ meta, old_excel, new_excel -> [ old_excel, new_excel ] }, 
                    griphins_ch.map{ meta, old_excel, new_excel -> meta.full_project_id },
                    [],
                    params.coverage,
                    params.bldb,
                    true,
                    griphins_ch.map{ meta, old_excel, new_excel -> meta.project_id }
                )
                ch_versions = ch_versions.mix(UPDATE_GRIPHIN.out.versions)

                // RESTORED: Assigning output variables for downstream use
                griphin_tsv_report = UPDATE_GRIPHIN.out.griphin_tsv_report
                griphin_report     = UPDATE_GRIPHIN.out.griphin_report

                // Example for one branch:
                // Join the versions with our 'collected' channel to get the 'thick' meta
                software_versions_ch = griphins_ch
                    .combine(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
                    .map { it ->
                        def meta_obj = it[0]
                        def v_file = it[-1]
                        // FIX: Include full_project_id from the original meta_obj
                        return [ [
                            id: meta_obj.project_id, 
                            project_id: meta_obj.project_id, 
                            full_project_id: meta_obj.full_project_id // Preserve the path we fixed earlier!
                        ], v_file ]
                    }

                META_CUSTOM_DUMPSOFTWAREVERSIONS ( software_versions_ch )
            } else { 
                // BRANCH 2: --input WITH --outdir (GRIPHIN_PUBLISH aggregated)
                // Extract old versions from collected_summaries_ch (which has the info files)
                ch_old_versions = collected_summaries_ch
                    .map { meta, files, path, busco, info -> 
                        def version = "unknown"
                        if (info instanceof Path && info.exists()) {
                            def content = info.text
                            if (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/) {
                                version = (content =~ /(?m)cdcgov\/phoenix:\s*(.+)/)[0][1].trim()
                            }
                        }
                        return version
                    }
                    .collect()
                    .map { versions -> 
                        def unique_versions = versions.unique()
                        if (unique_versions.size() > 1) {
                            println "WARNING (input with outdir): Multiple Phoenix versions detected: ${unique_versions}"
                        }
                        unique_versions.join(',')
                    }

                ch_inferred_mode = CREATE_INPUT_CHANNELS.out.mode_type.first()
                    .map { row ->
                        def (meta, mode_type) = row
                        mode_type.base == 'CDC_PHOENIX' || mode_type.base == 'PHOENIX' ? "READS" : "SCAFFOLD_INFERRED"
                    }

                ch_busco_check = collected_summaries_ch.map{ it[3] }.collect().map{ it.any{ it == true } }.ifEmpty(false)
                
                ch_clean_metas = griphin_inputs_ch
                    .filter { meta, file -> file != null && file.toString() != '[]' }
                    .combine(isolates_needing_update)
                    .filter { meta, file, update_ids -> 
                        update_ids != ['__EMPTY__'] && update_ids.contains(meta.id) 
                    }
                    .map { meta, file, update_ids -> [meta, file] }
                    .groupTuple(by: 0)
                    .map { m, f -> m + [filenames: f.flatten().findAll { it != null }.unique().collect { it.name }] }
                    .collect().ifEmpty([])

                ch_clean_files = griphin_inputs_ch.map{ it[1] }.collect().map{ it.flatten().unique() }.ifEmpty([])
                
                ch_output_dir = Channel.fromPath(params.outdir, type: 'dir').collect()

                GRIPHIN_PUBLISH (
                    params.ardb,
                    CREATE_INPUT_CHANNELS.out.valid_samplesheet.collect().ifEmpty([]),
                    ch_clean_metas,
                    ch_clean_files,
                    ch_output_dir,
                    workflow.manifest.version,
                    params.coverage,
                    ch_busco_check,        
                    ch_shigapass_check, 
                    false, 
                    params.bldb, 
                    true, 
                    false, 
                    [], 
                    ch_old_versions,      // ← NOW EXTRACTED FROM FILES
                    ch_inferred_mode      // ← NOW DETERMINED FROM DATA
                )
                
                ch_versions = ch_versions.mix(GRIPHIN_PUBLISH.out.versions)
                griphin_tsv_report = GRIPHIN_PUBLISH.out.griphin_tsv_report
                griphin_report     = GRIPHIN_PUBLISH.out.griphin_report

                // 1. Create the single versions file
                ch_collated_versions = ch_versions
                    .unique()
                    .collectFile(name: 'collated_versions.yml', sort: true)

                software_versions_ch = CREATE_INPUT_CHANNELS.out.directory_ch
                    .combine(isolates_needing_update)
                    .filter { meta, dir, update_ids -> 
                        update_ids != ['__EMPTY__'] && update_ids.contains(meta.id) 
                    }
                    .map { meta, dir, update_ids -> 
                        def p_id = meta.project_id.toString().split('/')[-1].replace("]", "")
                        [ [id: p_id, project_id: p_id, full_project_id: dir] ]
                    }
                    .unique()
                    .combine(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
                    .map { meta, v_file -> [meta, v_file] }

                META_CUSTOM_DUMPSOFTWAREVERSIONS (software_versions_ch)
            }
        }

    emit:
        mlst             = DO_MLST.out.checked_MLSTs
        amrfinder_output = AMRFINDERPLUS_RUN.out.report
        gamma_ar         = GAMMA_AR.out.gamma
        phx_summary      = GATHER_SUMMARY_LINES.out.summary_report
        //output for phylophoenix
        griphin_tsv      = griphin_tsv_report
        griphin_excel    = griphin_report
        //dir_samplesheet  = UPDATE_GRIPHIN.out.converted_samplesheet
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
