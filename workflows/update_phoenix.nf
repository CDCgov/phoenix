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
include { GAMMA as GAMMA_AR                    } from '../modules/local/gamma'
include { GET_TAXA_FOR_AMRFINDER               } from '../modules/local/get_taxa_for_amrfinder'
include { AMRFINDERPLUS_RUN                    } from '../modules/local/run_amrfinder'
include { CREATE_SUMMARY_LINE                  } from '../modules/local/phoenix_summary_line'
include { CREATE_AND_UPDATE_README             } from '../modules/local/updater/create_and_update_readme'
include { FETCH_FAILED_SUMMARIES               } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES                 } from '../modules/local/phoenix_summary'
include { GRIPHIN as GRIPHIN_NO_PUBLISH        } from '../modules/local/griphin'
include { GRIPHIN as GRIPHIN_NO_PUBLISH_CDC    } from '../modules/local/griphin'
include { GRIPHIN as GRIPHIN_PUBLISH           } from '../modules/local/griphin'
include { GRIPHIN as GRIPHIN_PUBLISH_CDC       } from '../modules/local/griphin'
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

include { CUSTOM_DUMPSOFTWAREVERSIONS as UPDATER_CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

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
    meta.project_id = folder_path.toString().split('/')[-1]
    def array = [ meta, griphin_ch ]
    return array
}

def add_meta_outdir(path_txt, griphin_excel, griphin_tsv, samplesheet) {
    // Step 1: Extract sample_name from griphin file
    def griphinFile = new File(griphin_tsv.toString())
    // Get second line, split on tab, get first element
    def sample_name = griphinFile.readLines()[1].split('\t')[0]
    //print("Sample name extracted from GRiPHin file: ${sample_name}\n")
    // Step 2: Extract the line with matching sample_name from samplesheet TSV file
    def matchingLine = null

    new File(samplesheet.toString()).withReader { reader ->
        // Skip header
        reader.readLine()
        // Read through the file line by line
        String line
        while ((line = reader.readLine()) != null) {
            // Split the line by tabs and check the first column
            def columns = line.split(',')
            if (columns[0].trim() == sample_name.trim()) {
                matchingLine = line
                print("Matching line found in samplesheet: ${matchingLine}\n")
                def meta = [:] // create meta array
                // Step 3: Extract the project_id from the matching line
                def cleanline = matchingLine.split(',')[1].trim()
                meta.project_id = cleanline.substring(0, cleanline.lastIndexOf('/'))
                //print("Meta project_id: ${meta.project_id}")
                return [ meta, griphin_excel ]
            }
        }
    }
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

def get_taxa_and_project_ID(input_ch){ 
    def genus = ""
    def species = ""
    input_ch[1].eachLine { line ->
        if (line.startsWith("G:")) {
            genus = line.split(":")[1].trim().split('\t')[1]
        } else if (line.startsWith("s:")) {
            species = line.split(":")[1].trim().split('\t')[1]
        }
    }
    //def clean_project_id = in_meta.project_id.replaceAll(/^['"]/, '').replaceAll(/['"]$/, '')
    return [input_ch[0], "$genus", input_ch[2] ]
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
        //filter to remove [[id:, project_id:], null, []] that comes out as part of remainder = true. Then either add [] if no readme is present or keep channel as is
        // Extract just the project name from the path
        files_to_update_ch = CREATE_INPUT_CHANNELS.out.pipeline_info.map{ meta, file -> [meta.project_id, meta, file] }
                                .combine(CREATE_INPUT_CHANNELS.out.directory_ch.map{ meta, dir -> [ meta.project_id.split('/').last(), meta, dir] }, by: [0])
                                .map { project_name, meta1, pipeline_info, meta2, directory_ch -> [meta2, pipeline_info, directory_ch]}
                                .join(CREATE_INPUT_CHANNELS.out.readme, by: [[0][0],[0][0]], remainder: true)
                                .filter{ it -> it.size() == 4 }.map{ meta, pipeline_info, dir, readme -> readme == null ? [meta, dir, pipeline_info, []] : [meta, dir, pipeline_info, readme] }

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

        ////////////////////////////////////// RUN SHIGAPASS IF IT WASN'T BEFORE AND IS CORRECT TAXA //////////////////////////////////////

        // First, create a set of isolate IDs that already have shigapass files --> we will use this to filter and only run samples that don't have shigapass files and need them
        existing_shigapass_ids = CREATE_INPUT_CHANNELS.out.shigapass.map{ meta, shigapass_file -> [ meta.id ]}.ifEmpty(["none"]) //[[id:], []]

        scaffolds_and_taxa_ch = CREATE_INPUT_CHANNELS.out.taxonomy.map{it -> get_taxa(it)}.filter{it, meta, taxonomy -> it.contains("Escherichia") || it.contains("Shigella")}.map{get_taxa_output, meta, taxonomy -> [[id:meta.id, project_id:meta.project_id], taxonomy ]}
                    .join(CREATE_INPUT_CHANNELS.out.filtered_scaffolds.map{                       meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}, by: [[0][0],[0][1]])
                    .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.map{ meta, fairy_outcome ->
                            // Read the content of the fairy_outcome file to check for FAILED
                            def content = file(fairy_outcome.toString()).text
                            def passed = !content.contains("FAILED")
                            [[id:meta.id, project_id:meta.project_id], [fairy_outcome, passed]]}, by: [[0][0],[0][1]])
                            .filter{ meta, taxonomy, filtered_scaffolds, fairy_data -> 
                                fairy_data[1] == true // Keep only entries where passed is true (no FAILED found)
                            }.map{   meta, taxonomy, filtered_scaffolds, fairy_data -> [meta, taxonomy, filtered_scaffolds] }.combine(existing_shigapass_ids)
                            .filter{ meta, taxonomy, filtered_scaffolds, existing_shigapass_ids -> existing_shigapass_ids == 'none' || !existing_shigapass_ids.contains(meta.id)}
                            .map{    meta, taxonomy, filtered_scaffolds, existing_shigapass_ids -> [meta, taxonomy, filtered_scaffolds ] }// Add the filter to exclude isolates that already have shigapass files

        /*/ For isolates that are E. coli or Shigella we will double check the FastANI Taxa ID and correct if necessary
        scaffolds_and_taxa_ch = CREATE_INPUT_CHANNELS.out.taxonomy.map{it -> get_taxa(it)}.filter{it, meta, taxonomy -> it.contains("Escherichia") || it.contains("Shigella")}.map{get_taxa_output, meta, taxonomy -> [[id:meta.id, project_id:meta.project_id], taxonomy ]}
            .join(CREATE_INPUT_CHANNELS.out.filtered_scaffolds.map{                               meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}, by: [[0][0],[0][1]])
            .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.splitCsv(strip:true, by:5).map{         meta, fairy_outcome      -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])
            .filter { meta, taxonomy, filtered_scaffolds, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}.combine(existing_shigapass_ids)
            .filter { meta, taxonomy, filtered_scaffolds, fairy_outcome, existing_shigapass_ids -> !existing_shigapass_ids.contains(meta.id)} // Add the filter to exclude isolates that already have shigapass files
            .map{ meta, taxonomy, filtered_scaffolds, fairy_outcome, existing_shigapass_ids -> return [meta, taxonomy, filtered_scaffolds ] }*/

        // Get ID from ShigaPass
        SHIGAPASS (
            scaffolds_and_taxa_ch, params.shigapass_database
        )
        ch_versions = ch_versions.mix(SHIGAPASS.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        checking_taxa_ch = CREATE_INPUT_CHANNELS.out.ani_best_hit.map{meta, ani_best_hit -> [[id:meta.id, project_id:meta.project_id], ani_best_hit]} \
            .join(CREATE_INPUT_CHANNELS.out.ani.map{                  meta, ani          -> [[id:meta.id, project_id:meta.project_id], ani ]},     by: [[0][0],[0][1]])\
            .join(SHIGAPASS.out.summary.map{                          meta, summary      -> [[id:meta.id, project_id:meta.project_id], summary ]}, by: [[0][0],[0][1]])

        // check shigapass and correct fastani taxa if its wrong
        CHECK_SHIGAPASS_TAXA (
            checking_taxa_ch
        )
        ch_versions = ch_versions.mix(CHECK_SHIGAPASS_TAXA.out.versions)

        ///////////////////////////////////// RUNNING AR CALLING //////////////////////////////////////

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file - written differently as fairy files might be variable in the number of lines in the file
        filtered_scaffolds_ch = CREATE_INPUT_CHANNELS.out.filtered_scaffolds.map{ meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}
                                .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.map{ meta, fairy_outcome ->
                                    // Read the content of the fairy_outcome file to check for FAILED
                                    def content = file(fairy_outcome.toString()).text
                                    def passed = !content.contains("FAILED")
                                    [[id:meta.id, project_id:meta.project_id], [fairy_outcome, passed]]}, by: [[0][0],[0][1]])
                                .filter { meta, filtered_scaffolds, fairy_data -> 
                                    fairy_data[1] // Keep only entries where passed is true (no FAILED found)
                                    }.map{ meta, filtered_scaffolds, fairy_data -> [meta, filtered_scaffolds] }

        /*/combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        filtered_scaffolds_ch = CREATE_INPUT_CHANNELS.out.filtered_scaffolds.map{    meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}
            .join(CREATE_INPUT_CHANNELS.out.fairy_outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome  -> [[id:meta.id, project_id:meta.project_id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])
            .filter { meta, filtered_scaffolds, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}
            .map{ meta, filtered_scaffolds, fairy_outcome -> return [meta, filtered_scaffolds] }*/

        // Running gamma to identify AR genes in scaffolds
        GAMMA_AR (
            filtered_scaffolds_ch, params.ardb
        )
        ch_versions = ch_versions.mix(GAMMA_AR.out.versions)

        // combing fastp_trimd information with fairy check of reads to confirm there are reads after filtering
        trimd_reads_file_integrity_ch = CREATE_INPUT_CHANNELS.out.reads.join(CREATE_INPUT_CHANNELS.out.fairy_outcome.map{ meta, fairy_outcome ->
                                    // Read the content of the fairy_outcome file to check for FAILED
                                    def content = file(fairy_outcome.toString()).text
                                    def passed = !content.contains("FAILED")
                                    [[id:meta.id, project_id:meta.project_id], [fairy_outcome, passed]]}, by: [[0][0],[0][1]])
                                .filter { meta, reads, fairy_data -> 
                                    fairy_data[1] == true // Keep only entries where passed is true (no FAILED found)
                                }.map{ meta, reads, fairy_data -> return [meta, reads] }
                                .join(CREATE_INPUT_CHANNELS.out.phoenix_tsv_ch.map{meta, phoenix_tsv_ch -> [[id:meta.id, project_id:meta.project_id], meta.entry]}, by: [[0][0],[0][1]])

        /*trimd_reads_file_integrity_ch = CREATE_INPUT_CHANNELS.out.reads.join(CREATE_INPUT_CHANNELS.out.fairy_outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])
                .filter { meta, reads, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}
                .map{ meta, reads, fairy_outcome -> return [meta, reads ] }
                .join(CREATE_INPUT_CHANNELS.out.phoenix_tsv_ch.map{meta, phoenix_tsv_ch -> [[id:meta.id, project_id:meta.project_id], meta.entry]}, by: [[0][0],[0][1]])*/

        // now we will split the channel into its true (busco present) and false (busco wasn't run with this dataset) elements
        busco_boolean_1ch = trimd_reads_file_integrity_ch.branch{ 
                    buscoTrue: it[3] == true
                    buscoFalse: it[3] == false}

        // Idenitifying AR genes in trimmed reads - using only datasets that were previously run with CDC_PHOENIX entry
        SRST2_AR (
            busco_boolean_1ch.buscoTrue.map{meta, reads, busco_boolean -> [meta, reads]}, "gene", params.ardb
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

        // Prepare channels with source information
        ch_ani_input_tagged = CREATE_INPUT_CHANNELS.out.ani_best_hit.map{ meta, file -> [meta.id, [meta, file, 'input']] }
        ch_ani_shigapass_tagged = CHECK_SHIGAPASS_TAXA.out.ani_best_hit.map{ meta, file -> [meta.id, [meta, file, 'shigapass']] }

        // For samples that were run through Shigapass we need to compare that to the list of previous ani_best_hit files, then only keep the shigapass files as that has the most up-to-date info
        // Mix and group by ID
        ch_ani_combined = ch_ani_input_tagged.mix(ch_ani_shigapass_tagged).groupTuple().map { id, values ->
                // First check if we have both input and shigapass files
                def shigapass_data = values.find { it[2] == 'shigapass' }
                def input_data = values.find { it[2] == 'input' }
                if (shigapass_data) {
                    // If shigapass data exists, use it
                    def (meta, file, _) = shigapass_data
                    return [meta, file]
                } else if (input_data) {
                    // Only use input data if no shigapass data
                    def (meta, file, _) = input_data
                    return [meta, file]
                } else {
                    // This shouldn't happen 
                    return null
                }
            }.filter { it != null } // Now ch_ani_combined contains the combined channel with prioritized files

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
        .join(ch_ani_combined.map{                                     meta, ani_best_hit    -> [[id:meta.id, project_id:meta.project_id], ani_best_hit]},   by: [[0][0],[0][1]])

        // First, check if SHIGAPASS.out.summary is empty and create appropriate channel - for cases where isolate set has no e.coli or shigella
        shigapass_ch = SHIGAPASS.out.summary.mix(CREATE_INPUT_CHANNELS.out.shigapass)

        // Extract [id, project_id] pairs for all isolates in line_summary_ch channel for joining // If shigapass_files is null (no match), replace with empty list
        all_id_pairs_ch = line_summary_ch.join(shigapass_ch, by: [[0][0],[0][1]], remainder: true).filter{ it[1] != null }
                .map{ meta, fastp_total_qc,checked_MLSTs,gamma_hv,gamma,gamma_pf,quast_report,assembly_ratio,synopsis,taxonomy,k2_bh_summary,report,ani_best_hit, shigapass_file ->
                if (shigapass_file == null) { return [meta, fastp_total_qc,checked_MLSTs,gamma_hv,gamma,gamma_pf,quast_report,assembly_ratio,synopsis,taxonomy,k2_bh_summary,report,ani_best_hit, []]
                } else { return [meta, fastp_total_qc,checked_MLSTs,gamma_hv,gamma,gamma_pf,quast_report,assembly_ratio,synopsis,taxonomy,k2_bh_summary,report,ani_best_hit, shigapass_file]}}

        line_summary_ch = all_id_pairs_ch

        // Generate summary per sample
        CREATE_SUMMARY_LINE (
            line_summary_ch
        )
        ch_versions = ch_versions.mix(CREATE_SUMMARY_LINE.out.versions)

        // Extract the list of project folders from input_dir channel, combine with the actual path that is needed by GATHER_SUMMARY_LINES
        project_ids = CREATE_INPUT_CHANNELS.out.directory_ch.map { it[1] }.unique().collect()
            .map { dirs -> dirs.collectEntries { dir -> def dirName = dir.tokenize('/').last()
            return [dirName, dir]}}

        // Create an empty channel as a fallback for when SRST2_AR doesn't run
        dummy_gene_results = Channel.empty()
        // Mix the real gene results with the dummy channel. This way, if SRST2_AR.out.gene_results doesn't exist, the empty channel is used
        def gene_results_ch
        try { gene_results_ch = SRST2_AR.out.gene_results.map{meta, gene_results -> [[id:meta.id, project_id:meta.project_id], gene_results] }
        } catch (Exception e) { gene_results_ch = dummy_gene_results }

        // Group the line_summaries by their project id and add in the full path for the project dir the join is only to make sure SRST2 finished before going to the last step, we just combine and then kick it out
        // remainder: true is used to ensure that if there are no gene results, the channel is still created
        summaries_ch = CREATE_SUMMARY_LINE.out.line_summary
            .join(gene_results_ch.map{meta, gene_results -> [[id:meta.id, project_id:meta.project_id], gene_results]}, by: [[0][0],[0][1]], remainder: true)
            .map{meta, summaryline, gene_results -> [meta.project_id, summaryline]}.groupTuple(by: 0)
            .map { group -> 
                def meta = [:] // create meta array
                def (id, files) = group
                id2 = id.split('/')[-1] // Remove full path and cut just the project id for use
                meta.project_id = id.split('/')[-1] // Remove full path and cut just the project id for use
                def dirPath = project_ids.value[id2]  // Retrieve the matching directory path
                return [meta, dirPath, files]}

        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            summaries_ch.map{ project_id, dir, summary_lines -> project_id}, 
            summaries_ch.map{ project_id, dir, summary_lines -> summary_lines}, 
            summaries_ch.map{ project_id, dir, summary_lines -> dir}, true
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        // check if the incoming file was created with CDC_PHOENIX or PHOENIX
        buscoExists_ch = summaries_ch.map{ project_id, dir, summary_lines -> 
                def busco_boolean = summary_lines.head().text.contains('BUSCO')
                return [ project_id, summary_lines, dir, busco_boolean ]}

        // now we will split the channel into its true and false elements for if we need to run the CDC version of GRIPHIN or not. 
        busco_boolean_ch = buscoExists_ch.branch{
                    buscoTrue: it[3] == true 
                    buscoFalse: it[3] == false}

        // If the input directory is not null, we need to check if the input directory is the same as the output directory
        if (params.outdir != "${launchDir}/phx_output") {
            outdir_full_path = Channel.fromPath(params.outdir)
            //this means that files need to be directed to the --outdir so we need to update the dir in the channel 
            buscoTrue_updated = busco_boolean_ch.buscoTrue.combine(outdir_full_path).map{project_id, summary_line, dir, busco_boolean, outdir -> [project_id, summary_line, outdir, busco_boolean] }
            buscoFalse_updated = busco_boolean_ch.buscoFalse.combine(outdir_full_path).map{project_id, summary_line, dir, busco_boolean, outdir -> [project_id, summary_line, outdir, busco_boolean] }

             // Use the updated channels
            busco_boolean_ch = [buscoTrue: buscoTrue_updated, buscoFalse: buscoFalse_updated]
        }

        // Check to see if the any isolates are Clostridioides difficile - set centar_var to true if it is, otherwise false
        // This is used to double check params.centar to ensure that griphin parameters are set correctly
        // collect all taxa and one by one count the number of c diff. then collect and get the sum to compare to 0
        centar_var = CREATE_INPUT_CHANNELS.out.taxonomy.map{ it -> get_only_taxa(it) }.collect().flatten().count{ it -> it == "Clostridioides"}.collect().sum().map{ it -> it[0] > 0 }
        // Now we need to check if --centar was passed when the samples were run previously. // "ifEmpty()" branch executes if no files match 
        centar_var = CREATE_INPUT_CHANNELS.out.centar.filter{ meta, file -> file.toString().endsWith("_centar_output.tsv") }
                        .ifEmpty { Channel.of(false)}.map { file -> return true}.first()  // If we get here with a file, it means we found a match and Take just the first result (true or false)

        if (params.indir != null) { // If the input directory is not null, we need to check if the input directory is the same as the output directory
            //pull in species specific files - use function to get taxa name, collect all taxa and one by one count the number of e. coli or shigella. then collect and get the sum to compare to 0
            shigapass_var = CREATE_INPUT_CHANNELS.out.taxonomy.map{it -> get_only_taxa(it)}.collect().flatten().count{ it -> it.contains("Escherichia") || it.contains("Shigella")}
                            .collect().sum().map{ it -> it[0] > 0 }

            // for samples that were created with entry CDC_PHX
            GRIPHIN_PUBLISH_CDC (
                busco_boolean_ch.buscoTrue.map{project_id, summary_line, dir, busco_boolean -> summary_line}, \
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, \
                busco_boolean_ch.buscoTrue.map{project_id, summary_line, dir, busco_boolean -> [dir.toString(), []]}, workflow.manifest.version, \
                params.coverage, false, shigapass_var, false, params.bldb, false
            )
            ch_versions = ch_versions.mix(GRIPHIN_PUBLISH_CDC.out.versions)

            // for samples that were created with entry PHX
            GRIPHIN_PUBLISH (
                busco_boolean_ch.buscoFalse.map{project_id, summary_line, dir, busco_boolean -> summary_line}, \
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, \
                busco_boolean_ch.buscoFalse.map{project_id, summary_line, dir, busco_boolean -> [dir.toString(), []]}, workflow.manifest.version, \
                params.coverage, true, shigapass_var, false, params.bldb, false
            )
            ch_versions = ch_versions.mix(GRIPHIN_PUBLISH.out.versions)

            griphin_tsv_report = GRIPHIN_PUBLISH_CDC.out.griphin_tsv_report.mix(GRIPHIN_PUBLISH.out.griphin_tsv_report)
            griphin_report = GRIPHIN_PUBLISH_CDC.out.griphin_report.mix(GRIPHIN_PUBLISH.out.griphin_report)

        } else {

            // check for shigapass
            shigapass_var_ch = CREATE_INPUT_CHANNELS.out.taxonomy.map{ meta, tax -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", "")], tax, meta.project_id]}
                                .map{it -> get_taxa_and_project_ID(it)}.groupTuple().map{meta, shiga_var, project_id -> [ meta, shiga_var, project_id.unique() ]}
                                .map{ meta, genus_list, project_id ->  [meta, genus_list.any{ genus -> genus == "Escherichia" || genus == "Shigella" }, project_id] }

            //add in the samplesheet specific to each project dir
            busco_true_ch = busco_boolean_ch.buscoTrue.join(CREATE_INPUT_CHANNELS.out.samplesheet_meta_ch, by: [0]).join(shigapass_var_ch, by: [0])
            busco_false_ch = busco_boolean_ch.buscoFalse.join(CREATE_INPUT_CHANNELS.out.samplesheet_meta_ch, by: [0]).join(shigapass_var_ch, by: [0])

            // for samples that were created with entry CDC_PHX
            GRIPHIN_NO_PUBLISH_CDC (
                busco_true_ch.map{meta, summary_line, dir, busco_boolean, samplesheet, shigapass_var, project_id -> summary_line}, \
                busco_true_ch.map{meta, summary_line, dir, busco_boolean, samplesheet, shigapass_var, project_id -> samplesheet}, \
                params.ardb, \
                busco_true_ch.map{meta, summary_line, dir, busco_boolean, samplesheet, shigapass_var, project_id -> [dir.toString(), project_id]}, workflow.manifest.version, \
                params.coverage, false, \
                busco_true_ch.map{meta, summary_line, dir, busco_boolean, samplesheet, shigapass_var, project_id -> shigapass_var}, \
                false, params.bldb, false
            )
            ch_versions = ch_versions.mix(GRIPHIN_NO_PUBLISH_CDC.out.versions)

            // for samples that were created with entry PHX
            GRIPHIN_NO_PUBLISH (
                busco_false_ch.map{meta, summary_line, dir, busco_boolean, samplesheet, shigapass_var, project_id -> summary_line}, \
                busco_false_ch.map{meta, summary_line, dir, busco_boolean, samplesheet, shigapass_var, project_id -> samplesheet}, \
                params.ardb, \
                busco_false_ch.map{meta, summary_line, dir, busco_boolean, samplesheet, shigapass_var, project_id -> [dir.toString(), project_id]}, workflow.manifest.version, \
                params.coverage, true, \
                busco_false_ch.map{meta, summary_line, dir, busco_boolean, samplesheet, shigapass_var, project_id -> shigapass_var}, \
                false, params.bldb, true
            )
            ch_versions = ch_versions.mix(GRIPHIN_NO_PUBLISH.out.versions)

            //define var to be used globally
            def griphin_reports_ch

            if (params.outdir == "${launchDir}/phx_output") {
                // If the output directory is not the default, we need to update the path in the channel
                griphin_reports_ch = GRIPHIN_NO_PUBLISH.out.griphin_report.collect().ifEmpty([]).combine(GRIPHIN_NO_PUBLISH_CDC.out.griphin_report.collect().ifEmpty([])).flatten().collate(2)
                                    .map{path_txt, griphin_report -> add_meta(path_txt, griphin_report)}
            } else {
                GRIPHIN_NO_PUBLISH.out.griphins.collect().ifEmpty([]).combine(GRIPHIN_NO_PUBLISH_CDC.out.griphins.collect().ifEmpty([])).flatten().collate(3).combine(CREATE_INPUT_CHANNELS.out.valid_samplesheet).view()
                griphin_reports_ch = GRIPHIN_NO_PUBLISH.out.griphins.collect().ifEmpty([]).combine(GRIPHIN_NO_PUBLISH_CDC.out.griphins.collect().ifEmpty([])).flatten().collate(3)
                                    .combine(CREATE_INPUT_CHANNELS.out.valid_samplesheet).map{path_txt, griphin_excel, griphin_tsv, samplesheet -> add_meta_outdir(path_txt, griphin_excel, griphin_tsv, samplesheet)}
            }

            // join old and new griphins for combining
            griphins_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch.map{meta, griphin_excel_ch -> [[project_id:meta.project_id], griphin_excel_ch]}.unique()\
                .join(griphin_reports_ch.map{                            meta, griphin_report   -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", "")], griphin_report]}, by: [0])\
                .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{        meta, directory_ch     -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", "")], directory_ch]}.unique(), by: [0])


            // combine griphin files, the new one just created and the old one that was found in the project dir. 
            UPDATE_GRIPHIN (
                griphins_ch.map{ meta, old_excel, new_excel, directory -> [ old_excel, new_excel ] }, 
                griphins_ch.map{ meta, old_excel, new_excel, directory -> directory },
                griphins_ch.map{ meta, old_excel, new_excel, directory -> meta.project_id}, // Get the project_id from the directory
                //griphins_ch.map{ meta, excel, report, directory, samplesheet -> samplesheet },
                [],
                params.coverage,
                params.bldb,
                true
            )
            ch_versions = ch_versions.mix(UPDATE_GRIPHIN.out.versions)

            griphin_tsv_report = UPDATE_GRIPHIN.out.griphin_tsv_report
            griphin_report = UPDATE_GRIPHIN.out.griphin_report
        }

        /// need to figure out how to do this on a per directory 
        // Collecting the software versions
        UPDATER_CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

    emit:
        mlst             = DO_MLST.out.checked_MLSTs
        amrfinder_output = AMRFINDERPLUS_RUN.out.report
        gamma_ar         = GAMMA_AR.out.gamma
        phx_summary      = GATHER_SUMMARY_LINES.out.summary_report
        //output for phylophoenix
        griphin_tsv      = griphin_tsv_report
        griphin_excel    = griphin_report
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
