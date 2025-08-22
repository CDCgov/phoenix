//
// Subworkflow: for updater and species specific phoenix summary and griphin output
//

include { GATHER_SUMMARY_LINES as SPECIES_SPECIFIC_GATHER_SUMMARY_LINES } from '../../modules/local/phoenix_summary' // calling it centar so output can stay together. 


workflow CREATE_PHOENIX_SUMMARY {
    take:
        line_summary     // CREATE_INPUT_CHANNELS.out.line_summary --> channel: tuple (meta.project_id, ) path(line_summary)
        directory_ch     // CREATE_INPUT_CHANNELS.out.directory_ch

    main:
        ch_versions = Channel.empty() // Used to collect the software versions


        // get summary lines and directory information to make sure all samples for a particular project folder stay together. 
        summaries_ch = line_summary.map{       meta, line_summary -> [[project_id:meta.project_id], line_summary] }
                        .join(directory_ch.map{meta, dir          -> [[project_id:meta.project_id], dir]}, by: [0])
                        .map{ meta, summary_line, dir -> [ meta, summary_line, dir, summary_line.readLines().first().contains('BUSCO') ]}

        transformed_summaries_ch = summaries_ch.map { meta, summary_line, dir, busco_boolean -> 
            def new_meta = [
                project_id: meta.project_id.toString().split('/')[-1].replace("]", ""), 
                full_project_id: dir  // assuming you want 'dir' as full_project_id
            ]
            return [new_meta, summary_line, dir, busco_boolean]}

        // For cases where the user is running isolates from multiple directories and wants the summary files output to their original project_id folders. 
        multiple_directories_ch = directory_ch.map{meta, dir -> [dir]}.collect().map{ files -> files.unique().size() > 1 }
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
            collected_summaries_ch = branched_collected_summaries_ch.ungrouped.collect().map{ items ->
                                def grouped_items = items.collate(4)
                                def summary_lines = []
                                def dirs = []
                                def busco_booleans = []
                            grouped_items.each{ item ->
                                summary_lines.add(item[1])
                                dirs.add(item[2])
                                busco_booleans.add(item[3])}
                            return [[], summary_lines, dirs.unique(), busco_booleans.any{ it == true }]}
        } else { 
            // Group by project and collect files within each project
            collected_summaries_ch = branched_collected_summaries_ch.ungrouped.mix(branched_collected_summaries_ch.grouped)
                .groupTuple(by: 0)  // group by meta (project_id)
                .map{ meta, summary_lines, dirs, busco_booleans -> 
                    [meta, summary_lines, dirs[0], busco_booleans[0]]}  // take first dir/busco since they should be the same within a project
        }

        
        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            gathered_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> meta},
            gathered_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> summary_lines},
            gathered_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> full_project_id},
            gathered_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> busco_boolean},
            gathered_summaries_ch.map{ meta, summary_lines, full_project_id, busco_boolean, pipeline_info -> pipeline_info}.map { file -> (file.text =~ /cdcgov\/phoenix: (.+)/)[0][1].trim() } // Extract the version from the pipeline_info file
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        // Combining sample summaries into final report
        SPECIES_SPECIFIC_GATHER_SUMMARY_LINES (
            collected_summaries_ch.map{ meta, summary_line, full_project_id, busco_boolean -> meta},
            collected_summaries_ch.map{ meta, summary_line, full_project_id, busco_boolean -> summary_line},
            collected_summaries_ch.map{ meta, summary_line,full_project_iddir, busco_boolean -> dir}, 
            collected_summaries_ch.map{ meta, summary_line, dir, busco_boolean -> busco_boolean},
            workflow.manifest.version
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

    emit:
        spades_ch                   = spades_ch
        spades_outcome              = SPADES.out.spades_outcome
        summary_line                = CREATE_SUMMARY_LINE_FAILURE.out.line_summary
        versions                    = ch_versions // channel: [ versions.yml ]
}