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
include { RENAME_FASTA_HEADERS           } from '../modules/local/rename_fasta_headers'
include { BUSCO                          } from '../modules/local/busco'
include { GAMMA_S as GAMMA_PF            } from '../modules/local/gammas'
include { GAMMA as GAMMA_AR              } from '../modules/local/gamma'
include { GAMMA as GAMMA_HV              } from '../modules/local/gamma'
include { BBMAP_REFORMAT                 } from '../modules/local/contig_less500'
include { SCAFFOLD_COUNT_CHECK           } from '../modules/local/fairy_scaffold_count_check'
include { QUAST                          } from '../modules/local/quast'
include { MASH_DIST                      } from '../modules/local/mash_distance'
include { FASTANI                        } from '../modules/local/fastani'
include { DETERMINE_TOP_MASH_HITS        } from '../modules/local/determine_top_mash_hits'
include { FORMAT_ANI                     } from '../modules/local/format_ANI_best_hit'
include { DETERMINE_TAXA_ID              } from '../modules/local/determine_taxa_id'
include { SHIGAPASS                      } from '../modules/local/shigapass'
include { CHECK_SHIGAPASS_TAXA           } from '../modules/local/check_shigapass_taxa'
include { PROKKA                         } from '../modules/local/prokka'
include { GET_TAXA_FOR_AMRFINDER         } from '../modules/local/get_taxa_for_amrfinder'
include { AMRFINDERPLUS_RUN              } from '../modules/local/run_amrfinder'
include { CALCULATE_ASSEMBLY_RATIO       } from '../modules/local/assembly_ratio'
include { CREATE_SUMMARY_LINE            } from '../modules/local/phoenix_summary_line'
include { FETCH_FAILED_SUMMARIES         } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES           } from '../modules/local/phoenix_summary'
include { GRIPHIN                        } from '../modules/local/griphin'

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/

include { CREATE_SCAFFOLDS_INPUT_CHANNEL           } from '../subworkflows/local/create_scaffolds_input_channel'
include { GENERATE_PIPELINE_STATS_WF               } from '../subworkflows/local/generate_pipeline_stats'
include { KRAKEN2_WF as KRAKEN2_ASMBLD             } from '../subworkflows/local/kraken2krona'
include { KRAKEN2_WF as KRAKEN2_WTASMBLD           } from '../subworkflows/local/kraken2krona'
include { DO_MLST                                  } from '../subworkflows/local/do_mlst'
include { CENTAR_SUBWORKFLOW                       } from '../subworkflows/local/centar_steps'

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

def add_empty_ch(input_ch) {
    meta_id = input_ch[0]
    output_array = [ meta_id, input_ch[1], input_ch[2], []]
    return output_array
}

def add_empty_chs(input_ch, int number_empty_ch) {
    // get number of items in the input ch
    int ch_size = input_ch.size() - 1 //subtract one to not take into account meta information
    meta_id = input_ch[0]
    // empty channel to add to
    def grouped_inputs = []
    grouped_inputs.add(meta_id)
    //loop through input channels and add to list
    for(int i in 1..ch_size){ grouped_inputs.add(input_ch[i]) }
    // create list of empty lists determined by input number
    for(int i in 1..number_empty_ch){ grouped_inputs.add([]) }
    // get correct channel size
    int complete_ch_size = ch_size + number_empty_ch
    output_array = grouped_inputs[0..complete_ch_size]
    return output_array
}

// Groovy funtion to make [ meta.id, [] ] - just an empty channel
def create_empty_ch(input_for_meta) { // We need meta.id associated with the empty list which is why .ifempty([]) won't work
    meta_id = input_for_meta[0]
    output_array = [ meta_id, [] ]
    return output_array
}

def get_taxa(input_ch){ 
        def genus = ""
        def species = ""
        input_ch[1].eachLine { line ->
            if (line.startsWith("G:")) {
                genus = line.split(":")[1].trim().split('\t')[1]
            } else if (line.startsWith("s:")) {
                species = line.split(":")[1].trim().split('\t')[1]
            }
        }
        //return ["$genus $species", input_ch[0], input_ch[1]]
        return ["$genus", input_ch[0], input_ch[1]]
}

def get_only_taxa(input_ch){ 
        def genus = ""
        def species = ""
        input_ch[1].eachLine { line ->
            if (line.startsWith("G:")) {
                genus = line.split(":")[1].trim().split('\t')[1]
            } else if (line.startsWith("s:")) {
                species = line.split(":")[1].trim().split('\t')[1]
            }
        }
        //return [ "$genus $species" ]
        return [ "$genus" ]
}

def add_project_id(old_meta, input_ch, outdir_path){
    def meta = [:] // create meta array
    meta.id = old_meta.id
    meta.project_id = outdir_path
    return [meta, input_ch]
}

def check_params_var(species_bol, species_param) {
    // species_bol -> was the species in question in the dataset?
    // species_param - >did the user pass the argument to run the species specific modules?
    if (species_bol == true && species_param == true){
        return true
    } else if (species_bol == true && species_param == false) {
        return false
    } else {
        return false
    }
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SCAFFOLDS_EXQC {
    take:
        ch_input
        ch_input_indir
        centar_param

    main:
        ch_versions = Channel.empty() // Used to collect the software versions
        // Allow outdir to be relative
        outdir_path = Channel.fromPath(params.outdir, relative: true)
        // Allow relative paths for krakendb argument
        kraken2_db_path  = Channel.fromPath(params.kraken2db, relative: true)

        CREATE_SCAFFOLDS_INPUT_CHANNEL (
            ch_input_indir, ch_input, ch_versions
        )
        ch_versions = ch_versions.mix(CREATE_SCAFFOLDS_INPUT_CHANNEL.out.versions)

        //unzip any zipped databases
        ASSET_CHECK (
            params.zipped_sketch, params.custom_mlstdb, kraken2_db_path
        )
        ch_versions = ch_versions.mix(ASSET_CHECK.out.versions)

        // Rename scaffold headers
        RENAME_FASTA_HEADERS (
            CREATE_SCAFFOLDS_INPUT_CHANNEL.out.scaffolds_ch
        )
        ch_versions = ch_versions.mix(RENAME_FASTA_HEADERS.out.versions)

        // Removing scaffolds <500bp
        BBMAP_REFORMAT (
            RENAME_FASTA_HEADERS.out.renamed_scaffolds
        )
        ch_versions = ch_versions.mix(BBMAP_REFORMAT.out.versions)

        // Combine bbmap log with the fairy outcome file
        scaffold_check_ch = BBMAP_REFORMAT.out.log.map{ it -> add_empty_chs(it, 6)}

        // Checking that there are still scaffolds left after filtering
        SCAFFOLD_COUNT_CHECK (
            scaffold_check_ch, true, params.coverage, params.nodes, params.names
        )
        ch_versions = ch_versions.mix(SCAFFOLD_COUNT_CHECK.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        filtered_scaffolds_ch = BBMAP_REFORMAT.out.filtered_scaffolds.map{    meta, filtered_scaffolds -> [[id:meta.id], filtered_scaffolds]}
            // The merging here is slightly different than with CDC_PHOENIX as the scaffolds outcome is meta.id meta.single_end for the scaffolds entry and only meta.id for CDC_PHOENIX
            .join(SCAFFOLD_COUNT_CHECK.out.outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome      -> [[id:meta.id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0])
            .filter { meta, filtered_scaffolds, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}
            .map{ meta, filtered_scaffolds, fairy_outcome -> return [meta, filtered_scaffolds] }

        // Running gamma to identify hypervirulence genes in scaffolds
        GAMMA_HV (
            filtered_scaffolds_ch, params.hvgamdb
        )
        ch_versions = ch_versions.mix(GAMMA_HV.out.versions)

        // Running gamma to identify AR genes in scaffolds
        GAMMA_AR (
            filtered_scaffolds_ch, params.ardb
        )
        ch_versions = ch_versions.mix(GAMMA_AR.out.versions)

        GAMMA_PF (
            filtered_scaffolds_ch, params.gamdbpf
        )
        ch_versions = ch_versions.mix(GAMMA_PF.out.versions)

        // Getting Assembly Stats
        QUAST (
            filtered_scaffolds_ch
        )
        ch_versions = ch_versions.mix(QUAST.out.versions)

        if (params.busco_db_path != null) {
            // Allow relative paths for krakendb argument
            busco_db_path = Channel.fromPath(params.busco_db_path, relative: true)
            // Add in krakendb into the fasta channel so each fasta has a krakendb to go with it.
            busco_ch = BBMAP_REFORMAT.out.filtered_scaffolds.combine(busco_db_path)\
            .join(SCAFFOLD_COUNT_CHECK.out.outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0])
            .filter { meta, filtered_scaffolds, busco_db_path, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}
            .map{ meta, filtered_scaffolds, busco_db_path, fairy_outcome -> return [meta, filtered_scaffolds, busco_db_path] }
        } else {
            // passing empty channel for busco db to align with expected inputs for the module
            busco_ch = BBMAP_REFORMAT.out.filtered_scaffolds.map{ meta, scaffolds -> [ [id:meta.id, single_end:meta.single_end], scaffolds, []]}\
            .join(SCAFFOLD_COUNT_CHECK.out.outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome      -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0])
            .filter { meta, filtered_scaffolds, busco_db_path, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}
            .map{ meta, filtered_scaffolds, busco_db_path, fairy_outcome -> return [meta, filtered_scaffolds, busco_db_path] }
        }

        // Checking single copy genes for assembly completeness
        BUSCO (
            busco_ch, 'auto', []
        )
        ch_versions = ch_versions.mix(BUSCO.out.versions)

        // Checking for Contamination in assembly creating krona plots and best hit files
        KRAKEN2_ASMBLD (
            BBMAP_REFORMAT.out.filtered_scaffolds, SCAFFOLD_COUNT_CHECK.out.outcome,"asmbld", [], QUAST.out.report_tsv, ASSET_CHECK.out.kraken_db
        )
        ch_versions = ch_versions.mix(KRAKEN2_ASMBLD.out.versions)

        // Creating krona plots and best hit files for weighted assembly
        KRAKEN2_WTASMBLD (
            BBMAP_REFORMAT.out.filtered_scaffolds, SCAFFOLD_COUNT_CHECK.out.outcome,"wtasmbld", [], QUAST.out.report_tsv, ASSET_CHECK.out.kraken_db
        )
        ch_versions = ch_versions.mix(KRAKEN2_WTASMBLD.out.versions)

        // combine filtered scaffolds and mash_sketch so mash_sketch goes with each filtered_scaffolds file
        mash_dist_ch = filtered_scaffolds_ch.combine(ASSET_CHECK.out.mash_sketch)

        // Running Mash distance to get top 20 matches for fastANI to speed things up
        MASH_DIST (
            mash_dist_ch
        )
        ch_versions = ch_versions.mix(MASH_DIST.out.versions)

        // Combining mash dist with filtered scaffolds and the outcome of the scaffolds count check based on meta.id
        top_mash_hits_ch = MASH_DIST.out.dist.join(filtered_scaffolds_ch, by: [0])

        // Generate file with list of paths of top taxa for fastANI
        DETERMINE_TOP_MASH_HITS (
            top_mash_hits_ch
        )
        ch_versions = ch_versions.mix(DETERMINE_TOP_MASH_HITS.out.versions)

        // Combining filtered scaffolds with the top taxa list based on meta.id
        top_taxa_list_ch = BBMAP_REFORMAT.out.filtered_scaffolds.map{meta, reads           -> [[id:meta.id], reads]}\
        .join(DETERMINE_TOP_MASH_HITS.out.top_taxa_list.map{         meta, top_taxa_list   -> [[id:meta.id], top_taxa_list ]}, by: [0])\
        .join(DETERMINE_TOP_MASH_HITS.out.reference_dir.map{         meta, reference_dir   -> [[id:meta.id], reference_dir ]}, by: [0])

        // Getting species ID
        FASTANI (
            top_taxa_list_ch
        )
        ch_versions = ch_versions.mix(FASTANI.out.versions)

        // Reformat ANI headers
        FORMAT_ANI (
            FASTANI.out.ani
        )
        ch_versions = ch_versions.mix(FORMAT_ANI.out.versions)

        // Combining weighted kraken report with the FastANI hit based on meta.id
        best_hit_ch = KRAKEN2_WTASMBLD.out.k2_bh_summary.map{meta, k2_bh_summary          -> [[id:meta.id], k2_bh_summary]}\
            .join(FORMAT_ANI.out.ani_best_hit_to_check.map{  meta, ani_best_hit_to_check  -> [[id:meta.id], ani_best_hit_to_check ]}, by: [0]).map{ it -> add_empty_ch(it) }

        // Getting ID from either FastANI or if fails, from Kraken2
        DETERMINE_TAXA_ID (
            best_hit_ch, params.nodes, params.names
        )
        ch_versions = ch_versions.mix(DETERMINE_TAXA_ID.out.versions)

        ////////////////////////////////////// SHIGAPASS //////////////////////////////////////
        // For isolates that are E. coli or Shigella we will double check the FastANI Taxa ID and correct if necessary
        // For isolates that are E. coli or Shigella we will double check the FastANI Taxa ID and correct if necessary
        scaffolds_and_taxa_ch = DETERMINE_TAXA_ID.out.taxonomy.map{it -> get_taxa(it)}.filter{it, meta, taxonomy -> it.contains("Escherichia") || it.contains("Shigella")}.map{get_taxa_output, meta, taxonomy -> [[id:meta.id], taxonomy ]}
            .join(BBMAP_REFORMAT.out.filtered_scaffolds.map{                                  meta, filtered_scaffolds -> [[id:meta.id], filtered_scaffolds]}, by: [0])
            .join(SCAFFOLD_COUNT_CHECK.out.outcome.splitCsv(strip:true, by:5).map{            meta, fairy_outcome      -> [[id:meta.id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0])
            .filter { meta, taxonomy, filtered_scaffolds, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}
            .map{ meta, taxonomy, filtered_scaffolds, fairy_outcome -> return [meta, taxonomy, filtered_scaffolds ] }

        // Get ID from ShigaPass
        SHIGAPASS (
            scaffolds_and_taxa_ch, params.shigapass_database
        )
        ch_versions = ch_versions.mix(SHIGAPASS.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        checking_taxa_ch = FORMAT_ANI.out.ani_best_hit_to_check.map{meta, ani_best_hit_to_check -> [[id:meta.id], ani_best_hit_to_check]} \
            .join(FASTANI.out.ani.map{                              meta, ani                   -> [[id:meta.id], ani ]},     by: [0])\
            .join(SHIGAPASS.out.summary.map{                        meta, summary               -> [[id:meta.id], summary ]}, by: [0])

        // check shigapass and correct fastani taxa if its wrong
        CHECK_SHIGAPASS_TAXA (
            checking_taxa_ch
        )
        ch_versions = ch_versions.mix(CHECK_SHIGAPASS_TAXA.out.versions)

        ////////////////////////////////////// PHOENIX //////////////////////////////////////
        // Perform MLST steps on isolates (with srst2 on internal samples)
        // 3rd input would normally be FASTP_TRIMD.out.reads, but since we have no reads we just pass something to get the meta_id and create an empty channel
        // This is necessary as we need the meta.id information for the mid_srst2_ch in the subworkflow.
        DO_MLST (
            BBMAP_REFORMAT.out.filtered_scaffolds, \
            SCAFFOLD_COUNT_CHECK.out.outcome, \
            RENAME_FASTA_HEADERS.out.renamed_scaffolds.map{ it -> create_empty_ch(it) }, \
            DETERMINE_TAXA_ID.out.taxonomy, \
            ASSET_CHECK.out.mlst_db, \
            false, "original"  //no reads to run srst2 and is opposed to the "update" option.
        )
        ch_versions = ch_versions.mix(DO_MLST.out.versions)

        ////////////////////////////////////// CENTAR ////////////////////////////////////// -- waiting for completed validation for release in v2.3.0
        // Run centar if necessary

        //First, check if any isolates are Clostridioides difficile and filter those to go through the channel
        determine_taxa_ch = DETERMINE_TAXA_ID.out.taxonomy.map{it -> get_taxa(it)}.filter{it, meta, taxonomy -> it == "Clostridioides"}.map{get_taxa_output, meta, taxonomy -> [[id:meta.id], taxonomy ]}

        if (centar_param == true) { // don't run regardless of what the isolates if --centar isn't passed
            // centar subworkflow requires project_ID as part of the meta
            CENTAR_SUBWORKFLOW (
                DO_MLST.out.checked_MLSTs.combine(outdir_path).map{meta, mlst, outdir -> add_project_id(meta, mlst, outdir)},
                SCAFFOLD_COUNT_CHECK.out.outcome.combine(outdir_path).map{meta, fairy, outdir -> add_project_id(meta, fairy, outdir)},
                BBMAP_REFORMAT.out.filtered_scaffolds.combine(outdir_path).map{meta, scaffolds, outdir -> add_project_id(meta, scaffolds, outdir)},
                ASSET_CHECK.out.mlst_db,
                determine_taxa_ch.combine(outdir_path).map{meta, taxa, outdir -> add_project_id(meta, taxa, outdir)}
            )
            ch_versions = ch_versions.mix(CENTAR_SUBWORKFLOW.out.versions)
        }

        ////////////////////////////////////// PHOENIX //////////////////////////////////////
        // get gff and protein files for amrfinder+
        PROKKA (
            filtered_scaffolds_ch, [], []
        )
        ch_versions = ch_versions.mix(PROKKA.out.versions)

        /*/ Fetch AMRFinder Database
        AMRFINDERPLUS_UPDATE( )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)*/

        // Create file that has the organism name to pass to AMRFinder
        GET_TAXA_FOR_AMRFINDER (
            DETERMINE_TAXA_ID.out.taxonomy, false
        )
        ch_versions = ch_versions.mix(GET_TAXA_FOR_AMRFINDER.out.versions)

        // Combining taxa and scaffolds to run amrfinder and get the point mutations.
        amr_channel = BBMAP_REFORMAT.out.filtered_scaffolds.map{                 meta, reads          -> [[id:meta.id], reads]}\
            .join(GET_TAXA_FOR_AMRFINDER.out.amrfinder_taxa.splitCsv(strip:true).map{meta, amrfinder_taxa -> [[id:meta.id], amrfinder_taxa ]}, by: [0])\
            .join(PROKKA.out.faa.map{                                                meta, faa            -> [[id:meta.id], faa ]},            by: [0])\
            .join(PROKKA.out.gff.map{                                                meta, gff            -> [[id:meta.id], gff ]},            by: [0])

        // Run AMRFinder
        AMRFINDERPLUS_RUN (
            amr_channel, params.amrfinder_db
        )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

        // Combining determined taxa with the assembly stats based on meta.id
        assembly_ratios_ch = DETERMINE_TAXA_ID.out.taxonomy.map{meta, taxonomy   -> [[id:meta.id], taxonomy]}\
            .join(QUAST.out.report_tsv.map{                     meta, report_tsv -> [[id:meta.id], report_tsv]}, by: [0])

        // Calculating the assembly ratio and gather GC% stats
        CALCULATE_ASSEMBLY_RATIO (
            assembly_ratios_ch, params.ncbi_assembly_stats
        )
        ch_versions = ch_versions.mix(CALCULATE_ASSEMBLY_RATIO.out.versions)

        // gather all outputs from shigapass and format_ani to get the best hit for each sample - we will flatten it to go into the pipeline stats
        ani_best_hit_ch = CHECK_SHIGAPASS_TAXA.out.ani_best_hit.collect().concat(FORMAT_ANI.out.ani_best_hit.collect()).flatten().collate(2)

        GENERATE_PIPELINE_STATS_WF (
            [], \
            [], \
            [], \
            [], \
            [], \
            [], \
            RENAME_FASTA_HEADERS.out.renamed_scaffolds, \
            BBMAP_REFORMAT.out.filtered_scaffolds, \
            DO_MLST.out.checked_MLSTs, \
            GAMMA_HV.out.gamma, \
            GAMMA_AR.out.gamma, \
            GAMMA_PF.out.gamma, \
            QUAST.out.report_tsv, \
            BUSCO.out.short_summaries_specific_txt, \
            KRAKEN2_ASMBLD.out.report, \
            KRAKEN2_ASMBLD.out.krona_html, \
            KRAKEN2_ASMBLD.out.k2_bh_summary, \
            KRAKEN2_WTASMBLD.out.report, \
            KRAKEN2_WTASMBLD.out.krona_html, \
            KRAKEN2_WTASMBLD.out.k2_bh_summary, \
            DETERMINE_TAXA_ID.out.taxonomy, \
            ani_best_hit_ch, \
            CALCULATE_ASSEMBLY_RATIO.out.ratio, \
            AMRFINDERPLUS_RUN.out.mutation_report, \
            CALCULATE_ASSEMBLY_RATIO.out.gc_content, \
            true
        )
        ch_versions = ch_versions.mix(GENERATE_PIPELINE_STATS_WF.out.versions)

        // Creating empty channel that has the form [ meta.id, [] ] that can be passed as a blank below
        empty_ch = RENAME_FASTA_HEADERS.out.renamed_scaffolds.map{ it -> create_empty_ch(it) }

        // Combining output based on meta.id to create summary by sample -- is this verbose, ugly and annoying? yes, if anyone has a slicker way to do this we welcome the input.
        line_summary_ch = empty_ch.map{                                  meta, list            -> [[id:meta.id], list]}\
        .join(DO_MLST.out.checked_MLSTs.map{                             meta, checked_MLSTs   -> [[id:meta.id], checked_MLSTs]},  by: [0])\
        .join(GAMMA_HV.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},          by: [0])\
        .join(GAMMA_AR.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},          by: [0])\
        .join(GAMMA_PF.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},          by: [0])\
        .join(QUAST.out.report_tsv.map{                                  meta, report_tsv      -> [[id:meta.id], report_tsv]},     by: [0])\
        .join(CALCULATE_ASSEMBLY_RATIO.out.ratio.map{                    meta, ratio           -> [[id:meta.id], ratio]},          by: [0])\
        .join(GENERATE_PIPELINE_STATS_WF.out.pipeline_stats.map{         meta, pipeline_stats  -> [[id:meta.id], pipeline_stats]}, by: [0])\
        .join(DETERMINE_TAXA_ID.out.taxonomy.map{                        meta, taxonomy        -> [[id:meta.id], taxonomy]},       by: [0])\
        .join(empty_ch.map{                                              meta, list            -> [[id:meta.id], list]},           by: [0])\
        .join(AMRFINDERPLUS_RUN.out.report.map{                          meta, report          -> [[id:meta.id], report]},         by: [0])\
        .join(ani_best_hit_ch.map{                                       meta, ani_best_hit    -> [[id:meta.id], ani_best_hit]},   by: [0])

        // Create a combined channel that contains all IDs from both line_summary_ch and SHIGAPASS.out.summary and handle the case where SHIGAPASS.out.summary might be empty
        shigapass_combined_ch = filtered_scaffolds_ch.map{ meta, scaffolds -> [[id:meta.id], meta.id] }  // Transform to [[meta.id], meta.id] for joining
                    .join(SHIGAPASS.out.summary, by: 0, remainder: true)  // Join on first element (meta.id)
                    .map{ id, original_id, shigapass_file -> [id, shigapass_file ?: []]}  // If shigapass_file is null, use empty list

        // Combine actual SHIGAPASS entries with backup empty entries and join with the original line_summary_ch
        line_summary_ch = line_summary_ch.join(shigapass_combined_ch, by: [0]) 

        // Generate summary per sample that passed SPAdes
        CREATE_SUMMARY_LINE (
            line_summary_ch
        )
        ch_versions = ch_versions.mix(CREATE_SUMMARY_LINE.out.versions)

        // Collect all the summary files prior to fetch step to force the fetch process to wait
        summaries_ch = CREATE_SUMMARY_LINE.out.line_summary.map{ meta, line_summary -> [line_summary]}.collect()
        // collect the failed fairy summary lines
        fairy_summary_ch = SCAFFOLD_COUNT_CHECK.out.summary_line.collect().ifEmpty( [] )
        // if centar was run, pull in species specific files
        if (centar_param == true) { // don't run regardless of what the isolates if --centar isn't passed
            centar_files_ch = CENTAR_SUBWORKFLOW.out.consolidated_centar.map{ meta, consolidated_file -> consolidated_file}.collect().ifEmpty([])
        }
        // if shigapass was run, pull in species specific files
        shigapass_files_ch = SHIGAPASS.out.summary.map{ meta, summary -> summary}.collect().ifEmpty([])
        // pulling it all together
        if (centar_param == true) { // don't run regardless of what the isolates if --centar isn't passed
            all_summaries_ch = spades_failure_summaries_ch.combine(failed_summaries_ch).combine(summaries_ch).combine(fairy_summary_ch).combine(centar_files_ch).combine(shigapass_files_ch)
        } else {
            all_summaries_ch = spades_failure_summaries_ch.combine(failed_summaries_ch).combine(summaries_ch).combine(fairy_summary_ch).combine(shigapass_files_ch)
        }

        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            all_summaries_ch, outdir_path, true
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        // Check to see if the any isolates are Clostridioides difficile - set centar_var to true if it is, otherwise false
        // This is used to double check params.centar to ensure that griphin parameters are set correctly
        //collect all taxa and one by one count the number of c diff. then collect and get the sum to compare to 0
        centar_boolean = DETERMINE_TAXA_ID.out.taxonomy.map{ it -> get_only_taxa(it) }.collect().flatten().count{ it -> it == "Clostridioides"}.collect().sum().map{ it -> it[0] > 0 }
        // Now we need to check if --centar was passed, In this case it is centar entry and therefore would be true
        centar_var = centar_boolean.map{ it -> check_params_var(it, centar_param)}
        //pull in species specific files - use function to get taxa name, collect all taxa and one by one count the number of e. coli or shigella. then collect and get the sum to compare to 0
        shigapass_var = DETERMINE_TAXA_ID.out.taxonomy.map{it -> get_only_taxa(it)}.collect().flatten().count{ it -> it.contains("Escherichia") || it.contains("Shigella")}
            .collect().sum().map{ it -> it[0] > 0 }

        //create GRiPHin report
        GRIPHIN (
            all_summaries_ch, CREATE_SCAFFOLDS_INPUT_CHANNEL.out.valid_samplesheet, params.ardb, outdir_path, workflow.manifest.version, params.coverage, false, true, false, shigapass_var, centar_var, params.bldb, true
        )
        ch_versions = ch_versions.mix(GRIPHIN.out.versions)

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
        // add in phx output for multiqc modulea
        ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.report_tsv.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_WTASMBLD.out.report.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(BUSCO.out.short_summaries_specific_txt.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    emit:
        scaffolds        = BBMAP_REFORMAT.out.filtered_scaffolds
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
