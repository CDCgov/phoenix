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

def create_meta(file_path, input_ch) {
    def meta = [:] // create meta array
    meta.project_id = file_path.splitText().first().toString().trim()
    array = [ meta, input_ch ]
    return array
}

def check_isolate_count(folderPath, meta){
    // Get all directories in the specified folder
    def excludeDirs = [ "pipeline_info","centar_pipeline_info", "multiqc"]
    def directories = []
    new File(folderPath).eachDir { dir ->
        // Only add if not in the exclude list
        if (!excludeDirs.contains(dir.name)) {
            directories.add(dir.name)
        }
    }
    // Keep track of unmatched directories
    def unmatchedDirs = []
    // Check which directories don't match any sample name
    directories.each { dirName ->
        if (!sampleNames.contains(dirName)) {
            unmatchedDirs.add(dirName)
        }
    }
    // Print unmatched directory names if any
    if (unmatchedDirs.size() > 0) {
        println "WARNING: The following directories do not match any sample name:"
        unmatchedDirs.each { dirName ->
            println "  - ${dirName}"
        }
    }
    // Return true if some directories matched and some didn't. (meaning there's at least one match AND at least one non-match)
    return unmatchedDirs.size() > 0 && unmatchedDirs.size() < directories.size()
}

def determine_entry(meta, tsv){
    // default is phoenix entry
    def noBusco = true
    try {
        // Open and read the TSV file
        def tsvFile = tsv.toFile()
        if (tsvFile.exists()) {
            // Read the first line (header)
            def headerLine = tsvFile.withReader { reader -> reader.readLine() }
            // Split the header by tab to get column names
            def columnNames = headerLine.split("\t")
            // Check if any column name contains "BUSCO" return true.
            noBusco = columnNames.any { colName -> 
                !colName.toUpperCase().contains("BUSCO") 
            }
        }
    } catch (Exception e) {
        // Handle any exceptions that might occur
        println "Error processing file $tsv: ${e.message}. see determine_entry() function in centar entry."
    }
    return [[project_id:meta.project_id], noBusco]
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
    def project_id = input_ch[0]
    def genera_list = input_ch[1]
    // Check if Escherichia is in the list
    if (!genera_list.contains("Clostridioides")) {
        errors.add("Project ${project_id} does not contain Clostridioides. Found genera: ${genera_list.join(', ')}")
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

//        ch_input_indir.view { "ch_input_indir: $it" }
//        ch_input.view { "input: $it" }
//        println ch_input_indir.getClass()
//        error "Stopping execution: custom reason here"

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
        summaries_ch = CREATE_INPUT_CHANNELS.out.line_summary.map{meta, line_summary -> [[project_id:meta.project_id], line_summary] }
            .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{     meta, dir          -> [[project_id:meta.project_id], dir]}, by: [0])
            .map { it -> filter_out_meta(it) }

        // Combining sample summaries into final report
        CENTAR_GATHER_SUMMARY_LINES (
            summaries_ch.map{ summary_line, dir -> summary_line}, summaries_ch.map{ summary_line, dir -> dir}, true
        )
        ch_versions = ch_versions.mix(CENTAR_GATHER_SUMMARY_LINES.out.versions)

        // collect centar output and summary lines to make the griphin report, note that groupTuple will force waiting until all samples are complete before proceeding
        griphin_input_ch = CREATE_INPUT_CHANNELS.out.line_summary.map{meta, line_summary -> [[project_id:meta.project_id], line_summary] }.groupTuple(by: [0])\
            .join(CENTAR_SUBWORKFLOW.out.consolidated_centar.map{     meta, centar_file  -> [[project_id:meta.project_id], centar_file]}.groupTuple(by: [0]), by: [0])\
            .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{         meta, dir          -> [[project_id:meta.project_id], dir]}, by: [0])\
            .join(CREATE_INPUT_CHANNELS.out.griphin_tsv_ch.map{       meta, tsv          -> determine_entry(meta, tsv)}, by: [0])

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
                griphin_input_ch.map{meta, summary_line, centar_files, dir, busco_var -> [dir]}, 
                workflow.manifest.version, params.coverage, 
                griphin_input_ch.map{meta, summary_line, centar_files, dir, busco_var -> busco_var}, 
                false, false, false, true, params.bldb, false
            )
            ch_versions = ch_versions.mix(CENTAR_GRIPHIN_INDIR.out.versions)

            griphin_report = CENTAR_GRIPHIN_INDIR.out.griphin_report

            // to be able to create software_versions.yml 
            software_versions_ch = ch_versions.unique().collectFile(name: 'collated_versions.yml')

            griphin_report.view { "griphin_report: $it" }

            // Combine with griphin_report to delay execution
            software_versions_ch
                .combine(griphin_report)
                .map { version, path_file, excel_file ->  // <-- Fixed: 3 parameters now
                    def project_dir = path_file.text.trim()  // <-- Read contents of path file
                    tuple(version, file(project_dir), project_dir.toString())
                }
                .set { software_versions_triplets }
            
            CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS(
                software_versions_triplets.map { v, p, s -> v },
                software_versions_triplets.map { v, p, s -> p },
                software_versions_triplets.map { v, p, s -> s }
            )

        } else if (params.input !=null) { // --input

            ///////////////////////////////// single dirs all samples in dir /////////////////////////////////////////////////
            // channel for isolates from the same dir, but we still need to check if we are running all the samples in the dir or not
            griphin_input_single_dir_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple.toBoolean() == false}
            filtering_samples_ch = griphin_input_single_dir_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> check_isolate_count(meta, dir) }
            griphin_input_single_all_isolates_ch = griphin_input_single_dir_ch.combine(filtering_samples_ch).map{meta, summary_lines, centar_files, dir, busco_var, is_multiple, filtering_samples -> filtering_samples.toBoolean() == false}

            //Run this process if there is only a single dir in the --input samples (via griphin_input_single_dir_ch) 
            //if no --outdir is not passed then isolates should go back to their project_ID folders (coded in modules.conf) and we want all samples in samplesheet to be included in the griphin report.
            //run CENTAR and use --filter_samples in griphin module (last var in channel input so only things in input are run) to make sure only samples in the --input are run
            CENTAR_GRIPHIN_INPUT (
                griphin_input_single_all_isolates_ch.map{meta, summary_line, centar_files, dir, busco_var, is_multiple, filtering_samples -> [summary_line, centar_files].flatten()}, 
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb, 
                griphin_input_single_all_isolates_ch.map{meta, summary_line, centar_files, dir, busco_var, is_multiple, filtering_samples -> [dir]}, 
                workflow.manifest.version, params.coverage, 
                griphin_input_single_all_isolates_ch.map{meta, summary_line, centar_files, dir, busco_var, is_multiple, filtering_samples -> busco_var}, 
                false, false, false, true, params.bldb, false
            )
            ch_versions = ch_versions.mix(CENTAR_GRIPHIN_INPUT.out.versions)

            griphin_report = CENTAR_GRIPHIN_INPUT.out.griphin_report

            ///////////////////////////////// single dirs some samples in dir /////////////////////////////////////////////////

            griphin_input_single_some_isolates_ch = griphin_input_single_dir_ch.combine(filtering_samples_ch).map{meta, summary_lines, centar_files, dir, busco_var, is_multiple, filtering_samples -> filtering_samples.toBoolean() == true}
            //griphin_input_single_some_isolates_ch.view()

            //Run this process if there is only a single dir in the --input samples (via griphin_input_single_some_isolates_ch)
            // if multiple_directories is true, then --input must have been used with samples from different directories. 
            // Also, if no --outdir is not passed then isolates should go back to their project_ID folders and we want all samples in that dir to be included in the griphin report. thus, we first make the griphin file with only the samples in --input
            NO_PUB_CENTAR_GRIPHIN (
                griphin_input_single_some_isolates_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple, filtering_samples -> [summary_lines, centar_files].flatten()}, 
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb,
                griphin_input_single_some_isolates_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple, filtering_samples -> [dir]}, 
                workflow.manifest.version, params.coverage,
                griphin_input_single_some_isolates_ch.map{meta, summary_line, centar_files, dir, busco_var, is_multiple, filtering_samples -> busco_var}, 
                false, true, false, true, params.bldb, true
            )
            ch_versions = ch_versions.mix(NO_PUB_CENTAR_GRIPHIN.out.versions)

            // collect centar output and summary lines to make the griphin report
            update_griphin_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch.map{meta, old_griphin_excel          -> [[project_id:meta.project_id], old_griphin_excel]}\
                .join(NO_PUB_CENTAR_GRIPHIN.out.griphin_report.map{            path_file, no_pub_griphin_report -> create_meta(path_file, no_pub_griphin_report)}, by: [0])\
                .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{              meta, dir                        -> [[project_id:meta.project_id], dir]}, by: [0])

            // combine original griphin file with the one that was just made, the new one just created and the old one that was found in the project dir. 
            UPDATE_CENTAR_GRIPHIN (
                update_griphin_ch.map{meta, old_griphin_excel, griphin_excel, dir -> [old_griphin_excel, griphin_excel]},
                update_griphin_ch.map{meta, old_griphin_excel, griphin_excel, dir -> dir},
                update_griphin_ch.map{meta, old_griphin_excel, griphin_excel, dir -> [dir]}.map{dir -> dir.toString().split('/')[-1].replace("]","")},
                CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                params.coverage,
                params.bldb,
                false
            )
            ch_versions = ch_versions.mix(UPDATE_CENTAR_GRIPHIN.out.versions)

            //griphin_report = UPDATE_CENTAR_GRIPHIN.out.griphin_report

            ///////////////////////////////// multiple dirs //////////////////////////////////////////////////////////////////
            // For cases where the user is running isolates from multiple directories and wants the summary files output to their original project_id folders. 
            multiple_directories_ch = CREATE_INPUT_CHANNELS.out.directory_ch.map{meta, dir -> [dir]}.collect().map{ files -> files.unique().size() > 1 }
            griphin_input_multi_dir_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple.toBoolean() == true}

            //Run this process if there is only a single dir in the --input samples (via griphin_input_multi_dir_ch)
            // if multiple_directories is true, then --input must have been used with samples from different directories. 
            // Also, if no --outdir is not passed then isolates should go back to their project_ID folders and we want all samples in that dir to be included in the griphin report. thus, we first make the griphin file with only the samples in --input
            NO_PUB_CENTAR_GRIPHIN_MULTI_DIR (
                griphin_input_multi_dir_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> [summary_lines, centar_files].flatten()}, 
                CREATE_INPUT_CHANNELS.out.valid_samplesheet, params.ardb,
                griphin_input_multi_dir_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> [dir]}, 
                workflow.manifest.version, params.coverage,
                griphin_input_multi_dir_ch.map{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> busco_var}, 
                false, true, false, true, params.bldb, true
            )
            ch_versions = ch_versions.mix(NO_PUB_CENTAR_GRIPHIN_MULTI_DIR.out.versions)

                        // collect centar output and summary lines to make the griphin report
            update_griphin_multi_dir_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch.map{meta, old_griphin_excel          -> [[project_id:meta.project_id], old_griphin_excel]}\
                .join(NO_PUB_CENTAR_GRIPHIN_MULTI_DIR.out.griphin_report.map{            path_file, no_pub_griphin_report -> create_meta(path_file, no_pub_griphin_report)}, by: [0])\
                .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{                        meta, dir                        -> [[project_id:meta.project_id], dir]}, by: [0])

            if (params.griphin_out == "${launchDir}" && params.outdir == "${launchDir}/phx_output") {

                // combine original griphin file with the one that was just made, the new one just created and the old one that was found in the project dir. 
                UPDATE_CENTAR_GRIPHIN_MULTI_DIR (
                    update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, dir -> [old_griphin_excel, griphin_excel]},
                    update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, dir -> dir},
                    update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, dir -> [dir]}.map{dir -> dir.toString().split('/')[-1].replace("]","")},
                    CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                    params.coverage,
                    params.bldb,
                    true
                )
                ch_versions = ch_versions.mix(UPDATE_CENTAR_GRIPHIN_MULTI_DIR.out.versions)

            } else {
                // Allow outdir to be relative
                griphin_out_path = Channel.fromPath(params.griphin_out.replace("./",""), relative: true)
                update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, dir -> [griphin_excel]}.flatten().collect()

                // combine original griphin file with the one that was just made, the new one just created and the old one that was found in the project dir. 
                UPDATE_CENTAR_GRIPHIN_MULTI_DIR (
                    update_griphin_multi_dir_ch.map{meta, old_griphin_excel, griphin_excel, dir -> [griphin_excel]}.flatten().collect(),
                    [],
                    griphin_out_path.map{dir -> dir.toString().split('/')[-1].replace("]","")},
                    CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                    params.coverage,
                    params.bldb,
                    false
                )
                ch_versions = ch_versions.mix(UPDATE_CENTAR_GRIPHIN_MULTI_DIR.out.versions)
            }

            griphin_report = UPDATE_CENTAR_GRIPHIN_MULTI_DIR.out.griphin_report.mix(UPDATE_CENTAR_GRIPHIN.out.griphin_report).mix(CENTAR_GRIPHIN_INPUT.out.griphin_report)

            if (params.outdir != "${launchDir}/phx_output") {
                // Single outdir: just emit to that one
                software_versions_ch = ch_versions.unique().collectFile(name: 'collated_versions.yml')

                software_versions_ch
                    .combine(outdir_path)
                    .map { version_file, dir_path ->
                        def path_obj = file(dir_path)
                        def val_obj = dir_path.toString()
                        tuple(version_file, path_obj, val_obj)
                    }
                    .set { software_versions_per_dir }

                CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS(
                    software_versions_per_dir.map { v, p, s -> v },
                    software_versions_per_dir.map { v, p, s -> p },
                    software_versions_per_dir.map { v, p, s -> s }
                )

            } else {
                // Multiple project output dirs
                software_versions_ch = ch_versions.unique().collectFile(name: 'collated_versions.yml')
                project_dirs_ch = update_griphin_ch
                    .map{meta, old_griphin_excel, griphin_excel, dir -> dir}
                    .mix(update_griphin_multi_dir_ch
                        .map{meta, old_griphin_excel, griphin_excel, dir -> dir})
                    .distinct()
                software_versions_ch
                    .combine(project_dirs_ch)
                    .map { version_file, dir_path ->
                        def path_obj = file(dir_path)
                        def val_obj = dir_path.toString()
                        tuple(version_file, path_obj, val_obj)
                    }
                    .set { software_versions_per_dir }
                CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS(
                    software_versions_per_dir.map { v, p, s -> v },
                    software_versions_per_dir.map { v, p, s -> p },
                    software_versions_per_dir.map { v, p, s -> s }
                )
            }
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