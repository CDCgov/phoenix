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

include { ASSET_CHECK                   } from '../modules/local/asset_check'
include { CREATE_SUMMARY_LINE           } from '../modules/local/phoenix_summary_line'
include { FETCH_FAILED_SUMMARIES        } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES          } from '../modules/local/phoenix_summary' // calling it centar so output can stay together. 
include { GRIPHIN as GRIPHIN_PUBLISH    } from '../modules/local/griphin' // calling it centar so output can stay together.
include { GRIPHIN as GRIPHIN_NO_PUBLISH } from '../modules/local/griphin'
include { UPDATE_GRIPHIN                } from '../modules/local/updater/update_griphin'

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

include { META_CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main_species_specific'

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

def add_meta(full_path_txt, griphin_ch) {
    def meta = [:] // create meta array
    def folder_path = full_path_txt.readLines().first()
    //meta.project_id = input_ch.getName().replaceAll("_GRiPHin.xlsx", "") // get file name without extention
    meta.project_id = folder_path.toString().split('/')[-1]
    def array = [ meta, griphin_ch ]
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
            params.zipped_sketch, params.custom_mlstdb, [], params.clia_amrfinder_db
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

        // For cases where the user is running isolates from multiple directories and wants the summary files output to their original project_id folders. 
        multiple_directories_ch = CREATE_INPUT_CHANNELS.out.directory_ch.map{meta, dir -> [dir]}.collect().map{ files -> files.unique().size() > 1 }
        // Conditionally group by project based on if there are multiple directories
        branched_collected_summaries_ch = transformed_summaries_ch
            .combine(multiple_directories_ch)
            .branch { meta, summary_line, dir, busco_boolean, is_multiple ->
                ungrouped: is_multiple == true
                    return [meta, summary_line, dir, busco_boolean]
                grouped: is_multiple == false
                    return [meta, summary_line, dir, busco_boolean] }

        if (params.outdir != "${launchDir}/phx_output") {
            //if there are multiple dirs and and --outdir given then the will need to all be combined into one Phoenix_Summary.tsv file in the outdir 
            collected_summaries_ungrouped_ch = branched_collected_summaries_ch.ungrouped.collect().map{ items ->
                                def grouped_items = items.collate(4)
                                def summary_lines = []
                                def dirs = []
                                def busco_booleans = []
                            grouped_items.each{ item ->
                                summary_lines.add(item[1])
                                dirs.add(item[2])
                                busco_booleans.add(item[3])}
                            return [[], summary_lines, dirs.unique(), busco_booleans.any{ it == true }]}
            //for single dirs
            collected_summaries_ch = branched_collected_summaries_ch.grouped.mix(collected_summaries_ungrouped_ch)
        } else { 
            // Group by project and collect files within each project
            collected_summaries_ch = branched_collected_summaries_ch.ungrouped.mix(branched_collected_summaries_ch.grouped)
                .groupTuple(by: 0)  // group by meta (project_id)
                .map{ meta, summary_lines, dirs, busco_booleans -> 
                    [meta, summary_lines, dirs[0], busco_booleans[0]]}  // take first dir/busco since they should be the same within a project
        }

        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> meta},
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> summary_line},
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> dir}, 
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> busco_boolean},
            workflow.manifest.version
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        // collect centar output and summary lines to make the griphin report, note that groupTuple will force waiting until all samples are complete before proceeding
        griphin_input_ch = CREATE_INPUT_CHANNELS.out.line_summary.map{meta, line_summary -> [[project_id:meta.project_id], line_summary] }.groupTuple(by: [0])\
            .join(CENTAR_SUBWORKFLOW.out.consolidated_centar.map{     meta, centar_file  -> [[project_id:meta.project_id], centar_file]}.groupTuple(by: [0]), by: [0])\
            .join(CREATE_INPUT_CHANNELS.out.griphin_tsv_ch.map{       meta, tsv          -> [[project_id:meta.project_id], tsv.readLines().first().contains('BUSCO')]}, by: [0])

        //define var to be used globally
        def griphin_report
        def griphin_out_path
        //griphin_input_multi_dir_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple.toBoolean() == true}
        // channel for isolates from the same dir, but we still need to check if we are running all the samples in the dir or not
        //griphin_input_single_dir_ch = griphin_input_ch.combine(multiple_directories_ch).filter{meta, summary_lines, centar_files, dir, busco_var, is_multiple -> is_multiple.toBoolean() == false}

        //create GRiPHin report channel
        griphin_inputs_ch = Channel.empty()
            .mix(
                CREATE_INPUT_CHANNELS.out.fastp_total_qc,
                CREATE_INPUT_CHANNELS.out.raw_stats,
                CREATE_INPUT_CHANNELS.out.k2_trimd_bh_summary,
                CREATE_INPUT_CHANNELS.out.k2_trimd_report,
                CREATE_INPUT_CHANNELS.out.k2_wtasmbld_bh_summary,
                CREATE_INPUT_CHANNELS.out.k2_wtasmbld_report,
                CREATE_INPUT_CHANNELS.out.quast_report,
                CREATE_INPUT_CHANNELS.out.fairy_outcome,
                CREATE_INPUT_CHANNELS.out.combined_mlst,
                CREATE_INPUT_CHANNELS.out.taxonomy,
                CREATE_INPUT_CHANNELS.out.assembly_ratio,
                CREATE_INPUT_CHANNELS.out.gc_content,
                CREATE_INPUT_CHANNELS.out.gamma_ar,
                CREATE_INPUT_CHANNELS.out.gamma_pf,
                CREATE_INPUT_CHANNELS.out.gamma_hv,
                CREATE_INPUT_CHANNELS.out.ani_best_hit,
                CREATE_INPUT_CHANNELS.out.synopsis,
                CREATE_INPUT_CHANNELS.out.busco_short_summary,
                CREATE_INPUT_CHANNELS.out.srst2_ar,
                CENTAR_SUBWORKFLOW.out.consolidated_centar
            )

        if (params.indir != null) { // --indir is passed then samples are from the same dir, thus all samples are run together and output in --indir unless --outdir or --griphin_out is passed

            busco_boolean = collected_summaries_ch.map{ meta, summary_lines, centar_file, busco_boolean -> busco_boolean}

            def outdir_full_path2
            if (params.outdir != "${launchDir}/phx_output") {
                outdir_full_path2 = Channel.fromPath(params.outdir, type: 'dir') // get the full path to the outdir, by not using "relative: true"
            } else {
                outdir_full_path2 = Channel.fromPath(params.indir, type: 'dir') // get the full path to the outdir, by not using "relative: true"
            }

            grouped_griphin_inputs_ch = griphin_inputs_ch.groupTuple()
                .map { meta, files ->
                    [
                        meta: [ id: "${meta.id}", filenames: files.collect { it.getName() } ],
                        files: files
                    ]
                }

            //create GRiPHin report
            GRIPHIN_PUBLISH (
                params.ardb,
                CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                grouped_griphin_inputs_ch.map { it.meta }.collect(),
                grouped_griphin_inputs_ch.map { it.files }.collect(),
                outdir_full_path2,
                workflow.manifest.version,
                params.coverage,
                // Fix once implemented, for now its hard-coded
                CREATE_INPUT_CHANNELS.out.griphin_tsv_ch.map{meta, tsv -> tsv.readLines().first().contains('BUSCO')}, 
                false, true, params.bldb, false, false, []// True is for centar_var but we are in centar so its always true" Falses are for filtering and 'dont_publish'
            )
            ch_versions = ch_versions.mix(GRIPHIN_PUBLISH.out.versions)

            griphin_tsv_report = GRIPHIN_PUBLISH.out.griphin_tsv_report
            griphin_report = GRIPHIN_PUBLISH.out.griphin_report

            // to be able to create software_versions.yml
            software_versions_ch = CREATE_INPUT_CHANNELS.out.directory_ch.map{meta, directory_ch -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", ""), full_project_id:directory_ch]]}.unique()
                                .combine(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

            META_CUSTOM_DUMPSOFTWAREVERSIONS(
                software_versions_ch
            )

        } else if (params.input !=null) { // --input

            if (params.outdir == "${launchDir}/phx_output") {

                //add in the samplesheet specific to each project dir
                boolean_ch = collected_summaries_ch.map{meta, summary_line, full_project_id, busco_boolean -> [[project_id:meta.project_id], full_project_id, busco_boolean]}
                                    .join(CREATE_INPUT_CHANNELS.out.samplesheet_meta_ch, by: [0])
                                    .map{meta, full_project_id, busco_boolean, samplesheet -> [[project_id:full_project_id], busco_boolean, samplesheet]}

                grouped_griphin_inputs_ch = griphin_inputs_ch.map { meta, file -> [meta.project_id, meta, file] }  // Restructure for grouping
                    // group all items by project_id
                    .groupTuple()
                    // build per-sample bundles inside each project
                    .map { String project_id, metas, files ->
                            def sample_ids = metas.unique{ it.id }*.id as List<String>
                            def per_sample = sample_ids.collect { sid ->
                                def sfiles = files.findAll { file ->
                                    def filename = file.getName()
                                    // Match files that start with sid followed by delimiter (explicit matching)
                                    filename.startsWith(sid + "_") || filename.startsWith(sid + ".") ||
                                    (filename.startsWith("short_summary") && filename.contains("." + sid + "."))
                                }
                                [
                                    meta : [ id: sid, filenames: sfiles.collect { it.getName() } ],
                                    files: sfiles
                                ]
                            }
                            def meta_list  = per_sample.collect { it.meta }
                            def files_flat = per_sample.collectMany { it.files }
                        // emit a tuple keyed by project for joining, plus the two payloads we need
                        [ [project_id: project_id], meta_list, files_flat ]
                    }

                // Join griphin_inputs_ch with boolean_ch
                combined_ch = grouped_griphin_inputs_ch.join(boolean_ch, by: [0])
                    .map { project_key, meta_list, files_flat, busco_boolean, samplesheet ->
                        // Order for the final call (we’ll unpack by index when invoking the module)
                        [
                            meta_list,              // 0 -> val(metas)
                            files_flat,             // 1 -> path(griphin_files)
                            project_key.project_id, // 2 -> path(outdir)
                            samplesheet,            // 3 -> path(original_samplesheet)
                            busco_boolean           // 4 -> val(entry)
                        ]
                    }

                //create GRiPHin report
                GRIPHIN_NO_PUBLISH (
                    params.ardb,                                   // path(db)
                    combined_ch.map { it[3] },                     // path(original_samplesheet)
                    combined_ch.map { it[0] },  // val(metas): list of [id:<sid>, filenames:[...]]
                    combined_ch.map { it[1] }, // path(griphin_files): flattened file list
                    combined_ch.map { it[2] },                     // path(outdir): full_project_id
                    workflow.manifest.version,                     // val(phx_version)
                    params.coverage,                               // val(coverage)
                    combined_ch.map { it[4] },                     // val(entry): busco_boolean
                    false,                                         // val(shigapass_detected)
                    true,                                          // val(centar_detected)
                    params.bldb,                                   // path(bldb)
                    true,                                          // val(filter_var)
                    true,                                           // val(dont_publish)
                    []
                )
                ch_versions = ch_versions.mix(GRIPHIN_NO_PUBLISH.out.versions)

                // If the output directory is not the default, we need to update the path in the channel
                griphin_reports_ch = GRIPHIN_NO_PUBLISH.out.griphin_report.collect().ifEmpty([[],[]]).flatten().collate(2)
                                    .map{path_txt, griphin_file -> add_meta(path_txt, griphin_file)}

                // join old and new griphins for combining
                griphins_ch = CREATE_INPUT_CHANNELS.out.griphin_excel_ch
                    .map{meta, griphin_excel_ch -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", "")], griphin_excel_ch]}.unique()
                    .join(griphin_reports_ch.map{                            meta, griphin_file   -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", "")], griphin_file]}, by: [0])
                    .join(CREATE_INPUT_CHANNELS.out.directory_ch.map{        meta, directory_ch   -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", "")], directory_ch]}.unique(), by: [0])
                    .map{meta, old_excel, new_excel, directory -> [[project_id:meta.project_id.toString().split('/')[-1].replace("]", ""), full_project_id:directory], old_excel, new_excel]}

                // combine griphin files, the new one just created and the old one that was found in the project dir. 
                UPDATE_GRIPHIN (
                    griphins_ch.map{ meta, old_excel, new_excel -> [ old_excel, new_excel ] }, 
                    griphins_ch.map{ meta, old_excel, new_excel -> meta.full_project_id },
                    //griphins_ch.map{ meta, excel, report, directory, samplesheet -> samplesheet },
                    [],
                    params.coverage,
                    params.bldb,
                    true,
                    griphins_ch.map{ meta, old_excel, new_excel -> meta.project_id }
                )
                ch_versions = ch_versions.mix(UPDATE_GRIPHIN.out.versions)

                griphin_tsv_report = UPDATE_GRIPHIN.out.griphin_tsv_report
                griphin_report = UPDATE_GRIPHIN.out.griphin_report

                software_versions_ch = griphins_ch.map{meta, old_excel, new_excel -> [[project_id:meta.project_id, full_project_id:meta.full_project_id]]}
                                .combine(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

                META_CUSTOM_DUMPSOFTWAREVERSIONS (
                    software_versions_ch
                )

            } else {

                outdir_full_path = Channel.fromPath(params.outdir, type: 'dir') // get the full path to the outdir, by not using "relative: true"
                //was busco run on any samples
                busco_boolean = summaries_ch.map{meta, summary_lines, full_project_id, busco_boolean -> busco_boolean}.collect().map{ busco_booleans -> busco_booleans.any{ it == true }}

                grouped_griphin_inputs_ch = griphin_inputs_ch.groupTuple()
                    .map { meta, files ->
                        [
                            meta: [ id: "${meta.id}", filenames: files.collect { it.getName() } ],
                            files: files 
                        ] }

                //create GRiPHin report
                GRIPHIN_PUBLISH (
                    params.ardb,
                    CREATE_INPUT_CHANNELS.out.valid_samplesheet,
                    grouped_griphin_inputs_ch.map { it.meta }.collect(),
                    grouped_griphin_inputs_ch.map { it.files }.collect(),
                    outdir_full_path,
                    workflow.manifest.version,
                    params.coverage,
                    busco_boolean,
                    false, true, params.bldb, true, false, []
                )
                ch_versions = ch_versions.mix(GRIPHIN_PUBLISH.out.versions)

                griphin_tsv_report = GRIPHIN_PUBLISH.out.griphin_tsv_report
                griphin_report = GRIPHIN_PUBLISH.out.griphin_report

                software_versions_ch = outdir_full_path.toList().map{dir -> [[project_id:dir.toString().split('/')[-1].replace("]", ""), full_project_id:dir]]}.unique()
                                .combine(ch_versions.unique().collectFile(name: 'collated_versions.yml'))

                META_CUSTOM_DUMPSOFTWAREVERSIONS (
                    software_versions_ch
                )
            }
        }


    emit:
        //output for phylophoenix
        griphins_excel   = griphin_report
        griphin_tsv      = griphin_tsv_report
        ch_versions      = software_versions_ch
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
