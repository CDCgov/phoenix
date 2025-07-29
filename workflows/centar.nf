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
include { GRIPHIN as CENTAR_GRIPHIN_INDIR                     } from '../modules/local/griphin' // calling it centar so output can stay together.
include { GRIPHIN as CENTAR_GRIPHIN_INPUT                     } from '../modules/local/griphin' // calling it centar so output can stay together.
include { GRIPHIN as NO_PUB_CENTAR_GRIPHIN                    } from '../modules/local/griphin'
include { GRIPHIN as NO_PUB_CENTAR_GRIPHIN_MULTI_DIR          } from '../modules/local/griphin' 
include { UPDATE_GRIPHIN as UPDATE_CENTAR_GRIPHIN             } from '../modules/local/updater/update_griphin'
include { UPDATE_GRIPHIN as UPDATE_CENTAR_GRIPHIN_MULTI_DIR   } from '../modules/local/updater/update_griphin'

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

include { CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main_centar'
include { CUSTOM_DUMPSOFTWAREVERSIONS        } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/


/*def filter_out_meta(input_ch) {
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
}*/

def create_meta(file_path, input_ch) {
    def meta = [:] // create meta array
    meta.project_id = file_path.splitText().first().toString().trim()
    array = [ meta, input_ch ]
    return array
}

def check_isolate_count(folderPath, centar_files){
    '''true means that not all samples in indir are being run through the pipeline.'''
    def excludeDirs = ["pipeline_info", "centar_pipeline_info", "multiqc"]
    def directories = file(folderPath.toString()).listFiles()
        ?.findAll { it.isDirectory() && !excludeDirs.contains(it.name) }
        ?.collect { it.name } ?: []
    // Extract sample IDs from centar and summary files
    def sampleIdsFromCentar = centar_files.collect{ file -> file.getName().replaceFirst(/_centar_output\.tsv$/, '') }
    def unmatchedDirs = directories.findAll { !sampleIdsFromCentar.contains(it) }
    //if (unmatchedDirs) { println "WARNING: Unmatched directories: ${unmatchedDirs.join(', ')}"}
    return unmatchedDirs.size() > 0 && unmatchedDirs.size() < directories.size()
}

def get_taxa_project_dir(input_ch){ 
        def genus = ""
        def species = ""
        input_ch[1].eachLine { line ->
            if (line.startsWith("G:")) {
                genus = line.split('\t')[1]
            } else if (line.startsWith("s:")) {
                species = line.split('\t')[1]
            }
        }
        return [input_ch[0], "$genus" ]
}

def validate_cdiff_presence(input_ch) {
    def errors = []
    // Process the flattened list in pairs (meta, genera_list)
    for (int i = 0; i < input_ch.size(); i += 2) {
        def project_meta = input_ch[i]
        def genera_list = input_ch[i + 1]
        //println("DEBUG: Processing project_meta: ${project_meta}")
        //println("DEBUG: Processing genera_list: ${genera_list}")
        
        // Extract project_id from the meta map
        def project_id = project_meta['project_id']
        // Convert ArrayBag to List first, then ensure all elements are strings
        genera_list = genera_list.toList().collect { it.toString() }
        // Check if Clostridioides is in the list
        if (!genera_list.contains("Clostridioides")) {
            errors.add("Project ${project_id} does not contain Clostridioides. Found genera: ${genera_list.join(', ')}")
        }
    }
    // If any errors found, exit with error message
    if (errors.size() > 0) {
        error("Pipeline validation failed:\n" + errors.join('\n'))
    }
    return true
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

        //confirm you have c diff samples in your run. 
        CREATE_INPUT_CHANNELS.out.taxonomy.map{meta, taxonomy -> [[project_id:meta.project_id], taxonomy] }.map{ it -> get_taxa_project_dir(it) }.groupTuple(by: [0])
                    .collect().map{ it -> validate_cdiff_presence(it)}

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
        summaries_ch = CREATE_INPUT_CHANNELS.out.line_summary.map{       meta, line_summary -> [[project_id:meta.project_id], line_summary] }
                        .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{meta, dir          -> [[project_id:meta.project_id], dir]}, by: [0])
                        .map{ meta, summary_line, dir -> [ meta, summary_line, dir, summary_line.readLines().first().contains('BUSCO') ]}

        transformed_summaries_ch = summaries_ch.map { meta, summary_line, dir, busco_boolean -> 
            def new_meta = [
                project_id: meta.project_id.toString().split('/')[-1].replace("]", ""), 
                full_project_id: dir  // assuming you want 'dir' as full_project_id
            ]
            return [new_meta, summary_line, dir, busco_boolean]
        }

        // Group by project and collect files within each project
        collected_summaries_ch = transformed_summaries_ch
            .groupTuple(by: 0)  // group by meta (project_id)
            .map{ meta, summary_lines, dirs, busco_booleans -> 
                [meta, summary_lines, dirs[0], busco_booleans[0]]  // take first dir/busco since they should be the same within a project
            }       

        collected_summaries_ch.view { "collated_summaries_ch: $it" }

        // Combining sample summaries into final report
        CENTAR_GATHER_SUMMARY_LINES (
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> meta},
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> summary_line},
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> dir}, 
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> busco_boolean}
        )
        ch_versions = ch_versions.mix(CENTAR_GATHER_SUMMARY_LINES.out.versions)

        

        // collect centar output and summary lines to make the griphin report, note that groupTuple will force waiting until all samples are complete before proceeding
        griphin_input_ch = CREATE_INPUT_CHANNELS.out.line_summary.map{meta, line_summary -> [[project_id:meta.project_id], line_summary] }.groupTuple(by: [0])\
            .join(CENTAR_SUBWORKFLOW.out.consolidated_centar.map{     meta, centar_file  -> [[project_id:meta.project_id], centar_file]}.groupTuple(by: [0]), by: [0])\
            .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{         meta, dir          -> [[project_id:meta.project_id], dir]}, by: [0])\
            .join(CREATE_INPUT_CHANNELS.out.griphin_tsv_ch.map{       meta, tsv          -> [[project_id:meta.project_id], !tsv.readLines().first().contains('BUSCO')]}, by: [0])


        //CREATE_INPUT_CHANNELS.out.griphin_tsv_ch.view()
        //[[project_id:/scicomp/groups-pure/OID/NCEZID/DHQP/CEMB/NCBI_HAISeq/Phoenix/haiseq_20250425], /scicomp/scratch/qpk9/3c/e3b17d1a66aa86c81360f6357325b8/haiseq_20250425_old_GRiPHin.tsv]
        /*griphin_input_ch.view() // for debugging purposes
        [[project_id:/scicomp/groups-pure/OID/NCEZID/DHQP/CEMB/NCBI_HAISeq/Phoenix/haiseq_20250425], 
        [/scicomp/scratch/qpk9/f7/d9b37b73d2b3082ee4fc7fbf6eb860/2025EL-00073/2025EL-00073_summaryline.tsv, 
        /scicomp/scratch/qpk9/2e/0fa99a6367fda163ccf47c3914c4ad/2025EL-00074/2025EL-00074_summaryline.tsv, 
        /scicomp/scratch/qpk9/04/776b0390373bca119c9bd63487aacf/2025EL-00075/2025EL-00075_summaryline.tsv, 
        /scicomp/scratch/qpk9/c3/1d41896210088b9e7880a02e9848dd/2025EL-00072/2025EL-00072_summaryline.tsv], 
        [/scicomp/scratch/qpk9/78/b75991fae1dc482595ed370a680286/2025EL-00074_centar_output.tsv, 
        /scicomp/scratch/qpk9/97/9738f7596617b9c9d854f29d15b12a/2025EL-00075_centar_output.tsv, 
        /scicomp/scratch/qpk9/46/5da030c72a8c3f07d3005bf1e37101/2025EL-00072_centar_output.tsv, 
        /scicomp/scratch/qpk9/5f/048037c0e6670cb0b4c4e5d9793876/2025EL-00073_centar_output.tsv], /scicomp/groups-pure/OID/NCEZID/DHQP/CEMB/NCBI_HAISeq/Phoenix/haiseq_20250425, true]*/

        //define var to be used globally
        def griphin_report
        def griphin_out_path
        // For cases where the user is running isolates from multiple directories and wants the summary files output to their original project_id folders. 
        multiple_directories_ch = CREATE_INPUT_CHANNELS.out.directory_ch.map{meta, dir -> [dir]}.collect().map{ files -> files.unique().size() > 1 }
        //griphin_input_multi_dir_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple.toBoolean() == true}
        // channel for isolates from the same dir, but we still need to check if we are running all the samples in the dir or not
        //griphin_input_single_dir_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple.toBoolean() == false}
        //griphin_input_single_dir_ch.view()

        if (params.indir != null) { // --indir is passed then samples are from the same dir, thus all samples are run together and output in --indir unless --outdir or --griphin_out is passed

            //create GRiPHin report
            CENTAR_GRIPHIN_INDIR (
                griphin_input_ch.map{meta, summary_line, centar_files, dir, busco_var -> [summary_line, centar_files].flatten()}, 
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, 
                griphin_input_ch.map{meta, summary_line, centar_files, dir, busco_var -> [dir, []]}, 
                workflow.manifest.version, params.coverage, 
                griphin_input_ch.map{meta, summary_line, centar_files, dir, busco_var -> busco_var}, 
                false, true, params.bldb, false
            )
            ch_versions = ch_versions.mix(CENTAR_GRIPHIN_INDIR.out.versions)

            griphin_report = CENTAR_GRIPHIN_INDIR.out.griphin_report

            // to be able to create software_versions.yml 
            software_versions_ch = ch_versions.unique().collectFile(name: 'collated_versions.yml')

            griphin_report.view { "griphin_report: $it" }

            //// Combine with griphin_report to delay execution
            // software_versions_ch
            //     .combine(griphin_report)
            //     .map { version, path_file, excel_file ->  // <-- Fixed: 3 parameters now
            //         def project_dir = path_file.text.trim()  // <-- Read contents of path file
            //         tuple(version, file(project_dir), project_dir.toString())
            //     }
            //     .set { software_versions_triplets }
            
            // CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS(
            //     software_versions_triplets.map { v, p, s -> v },
            //     software_versions_triplets.map { v, p, s -> p },
            //     software_versions_triplets.map { v, p, s -> s }
            //)
            CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS(
                software_versions_ch
            )

        } else if (params.input !=null) { // --input

            ///////////////////////////////// single dirs all samples in dir /////////////////////////////////////////////////
            // channel for isolates from the same dir, but we still need to check if we are running all the samples in the dir or not
            griphin_input_single_dir_to_split_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple.toBoolean() == false}
            filtering_samples_ch = griphin_input_single_dir_to_split_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> check_isolate_count(dir, centar_files) }

            // now we will split the channel into its true and false elements for if we are running all samples in the dir or only some 
            griphin_input_single_dir_ch = griphin_input_single_dir_to_split_ch.combine(filtering_samples_ch).branch{
                        single_some_isolates: it[6] == true 
                        single_all_isolates: it[6] == false}

            //Run this process if there is only a single dir in the --input samples (via griphin_input_single_dir_ch) 
            //if no --outdir is not passed then isolates should go back to their project_ID folders (coded in modules.conf) and we want all samples in samplesheet to be included in the griphin report.
            //run CENTAR and use --filter_samples in griphin module (last var in channel input so only things in input are run) to make sure only samples in the --input are run
            CENTAR_GRIPHIN_INPUT (
                griphin_input_single_dir_ch.single_all_isolates.map{meta, summary_line, centar_files, dir, busco_var, is_multiple, filtering_samples -> [summary_line, centar_files].flatten()}, 
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, 
                griphin_input_single_dir_ch.single_all_isolates.map{meta, summary_line, centar_files, dir, busco_var, is_multiple, filtering_samples -> [dir, []]}, 
                workflow.manifest.version, params.coverage, 
                griphin_input_single_dir_ch.single_all_isolates.map{meta, summary_line, centar_files, dir, busco_var, is_multiple, filtering_samples -> busco_var}, 
                false, true, params.bldb, false
            )
            ch_versions = ch_versions.mix(CENTAR_GRIPHIN_INPUT.out.versions)
            griphin_report = CENTAR_GRIPHIN_INPUT.out.griphin_report

            def update_griph_out_dir
            if (params.outdir == "${launchDir}/phx_output" && params.griphin_out != "${launchDir}") {
                update_griph_out_path = params.griphin_out
            } else {
                update_griph_out_path = params.outdir
            }
            update_griph_out_dir = update_griph_out_path.split('/')[-1]

            ///////////////////////////////// single dirs some samples in dir /////////////////////////////////////////////////

            //Run this process if there is only a single dir in the --input samples (via griphin_input_single_some_isolates_ch)
            // if multiple_directories is true, then --input must have been used with samples from different directories. 
            // Also, if no --outdir is not passed then isolates should go back to their project_ID folders and we want all samples in that dir to be included in the griphin report. thus, we first make the griphin file with only the samples in --input
            NO_PUB_CENTAR_GRIPHIN (
                griphin_input_single_dir_ch.single_some_isolates.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple, filtering_samples -> [summary_lines, centar_files].flatten()}, 
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb,
                griphin_input_single_dir_ch.single_some_isolates.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple, filtering_samples -> [dir, []]}, 
                workflow.manifest.version, params.coverage,
                griphin_input_single_dir_ch.single_some_isolates.map{meta, summary_line, centar_files, dir, busco_var, is_multiple, filtering_samples -> busco_var}, 
                false, true, params.bldb, true
            )
            ch_versions = ch_versions.mix(NO_PUB_CENTAR_GRIPHIN.out.versions)

            // collect centar output and summary lines to make the griphin report
            update_griphin_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch.map{meta, old_griphin_excel          -> [[project_id:meta.project_id], old_griphin_excel]}\
                .join(NO_PUB_CENTAR_GRIPHIN.out.griphin_report.map{            path_file, no_pub_griphin_report -> create_meta(path_file, no_pub_griphin_report)}, by: [0])\
                .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{              meta, dir                        -> [[project_id:meta.project_id], dir]}, by: [0])
                .map{meta, old_griphin_excel, no_pub_griphin_report, directory -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", ""), full_project_id:directory], old_excel, new_excel]}

            // combine original griphin file with the one that was just made, the new one just created and the old one that was found in the project dir. 
            UPDATE_CENTAR_GRIPHIN (
                update_griphin_ch.map{meta, old_griphin_excel, griphin_excel, dir -> [old_griphin_excel, griphin_excel]},
                update_griphin_ch.map{ meta, old_excel, new_excel, directory -> meta },
                CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                params.coverage,
                params.bldb,
                false,
                update_griph_out_dir
            )
            ch_versions = ch_versions.mix(UPDATE_CENTAR_GRIPHIN.out.versions)

            //griphin_report = UPDATE_CENTAR_GRIPHIN.out.griphin_report

            ///////////////////////////////// multiple dirs //////////////////////////////////////////////////////////////////
            // For cases where the user is running isolates from multiple directories and wants the summary files output to their original project_id folders. 
            //multiple_directories_ch = CREATE_INPUT_CHANNELS.out.directory_ch.map{meta, dir -> [dir]}.collect().map{ files -> files.unique().size() > 1 }
            multiple_directories_ch = CREATE_INPUT_CHANNELS.out.directory_ch
                .map{meta, dir -> dir}
                .collect()
                .map{ dirs -> dirs.unique().size() > 1 }
            //griphin_input_multi_dir_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple.toBoolean() == true}
            griphin_input_multi_dir_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple}

            //Run this process if there is only a single dir in the --input samples (via griphin_input_multi_dir_ch)
            // if multiple_directories is true, then --input must have been used with samples from different directories. 
            // Also, if no --outdir is not passed then isolates should go back to their project_ID folders and we want all samples in that dir to be included in the griphin report. thus, we first make the griphin file with only the samples in --input
            NO_PUB_CENTAR_GRIPHIN_MULTI_DIR (
                griphin_input_multi_dir_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> [summary_lines, centar_files].flatten()}, 
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb,
                griphin_input_multi_dir_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> [dir, []]}, 
                workflow.manifest.version, params.coverage,
                griphin_input_multi_dir_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> busco_var}, 
                false, true, params.bldb, true
            )
            ch_versions = ch_versions.mix(NO_PUB_CENTAR_GRIPHIN_MULTI_DIR.out.versions)

            // collect centar output and summary lines to make the griphin report
            update_griphin_multi_dir_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch.map{meta, old_griphin_excel          -> [[project_id:meta.project_id], old_griphin_excel]}\
                .join(NO_PUB_CENTAR_GRIPHIN_MULTI_DIR.out.griphin_report.map{            path_file, no_pub_griphin_report -> create_meta(path_file, no_pub_griphin_report)}, by: [0])\
                .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{                        meta, directory                        -> [[project_id:meta.project_id], directory]}, by: [0])
                .map{meta, old_griphin_excel, no_pub_griphin_report, directory -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", ""), full_project_id:directory], old_griphin_excel, no_pub_griphin_report, directory]}

            if (params.griphin_out == "${launchDir}" && params.outdir == "${launchDir}/phx_output") {

                // combine original griphin file with the one that was just made, the new one just created and the old one that was found in the project dir. 
                UPDATE_CENTAR_GRIPHIN_MULTI_DIR (
                    update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, directory -> [old_griphin_excel, griphin_excel]},
                    update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, directory -> meta },
                    CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                    params.coverage,
                    params.bldb,
                    true,
                    update_griph_out_dir
                )
                ch_versions = ch_versions.mix(UPDATE_CENTAR_GRIPHIN_MULTI_DIR.out.versions)

            } else {
                // Allow outdir to be relative
                griphin_out_path = Channel.fromPath(params.griphin_out.replace("./",""), relative: true)
                flattened_update_griphin_multi_dir_ch = update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, directory -> [griphin_excel]}.flatten().collect()
                flattened_meta = update_griphin_multi_dir_ch
                    .map{meta, old_griphin_excel, griphin_excel, directory -> meta}
                    .first()

                update_griphin_multi_dir_ch.view { "update_griphin_multi_dir_ch: $it" }



                // combine original griphin file with the one that was just made, the new one just created and the old one that was found in the project dir. 
                UPDATE_CENTAR_GRIPHIN_MULTI_DIR (
                    //update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, directory -> [griphin_excel]}.flatten().collect(),
                    flattened_update_griphin_multi_dir_ch,
                    //This will need to be updated, but right now it just grabs ONE of the meta's to pull info from.
                    flattened_meta,
                    CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                    params.coverage,
                    params.bldb,
                    false,
                    update_griph_out_dir
                )
                ch_versions = ch_versions.mix(UPDATE_CENTAR_GRIPHIN_MULTI_DIR.out.versions)
            }

            griphin_report = UPDATE_CENTAR_GRIPHIN_MULTI_DIR.out.griphin_report.mix(UPDATE_CENTAR_GRIPHIN.out.griphin_report).mix(CENTAR_GRIPHIN_INPUT.out.griphin_report)
            software_versions_ch = CREATE_INPUT_CHANNELS.out.directory_ch.map{meta, directory_ch -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", ""), full_project_id:directory_ch]]}.unique()
                                .combine(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

            // if (params.outdir != "${launchDir}/phx_output") {
            //     // Single outdir: just emit to that one
            //     software_versions_ch = ch_versions.unique().collectFile(name: 'collated_versions.yml')

            //     software_versions_ch
            //         .combine(outdir_path)
            //         .map { version_file, dir_path ->
            //             def path_obj = file(dir_path)
            //             def val_obj = dir_path.toString()
            //             tuple(version_file, path_obj, val_obj)
            //         }
            //         .set { software_versions_per_dir }

            //     CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS(
            //         software_versions_per_dir.map { v, p, s -> v },
            //         software_versions_per_dir.map { v, p, s -> p },
            //         software_versions_per_dir.map { v, p, s -> s }
            //     )

            // } else {
            //     // Multiple project output dirs
            //     software_versions_ch = ch_versions.unique().collectFile(name: 'collated_versions.yml')
            //     project_dirs_ch = update_griphin_ch
            //         .map{meta, old_griphin_excel, griphin_excel, dir -> dir}
            //         .mix(update_griphin_multi_dir_ch
            //             .map{meta, old_griphin_excel, griphin_excel, dir -> dir})
            //         .distinct()
            //     software_versions_ch
            //         .combine(project_dirs_ch)
            //         .map { version_file, dir_path ->
            //             def path_obj = file(dir_path)
            //             def val_obj = dir_path.toString()
            //             tuple(version_file, path_obj, val_obj)
            //         }
            //         .set { software_versions_per_dir }
                CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS(
                    software_versions_ch
//                    software_versions_per_dir.map { v, p, s -> v },
//                    software_versions_per_dir.map { v, p, s -> p },
//                    software_versions_per_dir.map { v, p, s -> s }
                )
//            }
        } else {
            exit 1, "You shouldn't be here, please open a github issue to report."
        }
        

    emit:
        //output for phylophoenix
        griphins_excel    = griphin_report
        //ch_versions       = software_versions_ch
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


/*
========================================================================================
    THE END
========================================================================================
*/