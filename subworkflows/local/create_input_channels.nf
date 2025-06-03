//
// workflow handles taking in either a samplesheet or directory and creates correct channels for scaffolds entry point
//
// for centar entry
include { SAMPLESHEET_CHECK as CENTAR_SAMPLESHEET_CHECK   } from '../../modules/local/samplesheet_check'
include { CREATE_SAMPLESHEET as CENTAR_CREATE_SAMPLESHEET } from '../../modules/local/create_samplesheet'
// for cdc_phoenix, phoenix entry
include { SAMPLESHEET_CHECK as SAMPLESHEET_CHECK          } from '../../modules/local/samplesheet_check'
include { CREATE_SAMPLESHEET as CREATE_SAMPLESHEET        } from '../../modules/local/create_samplesheet'
// for update entry
include { COLLECT_SAMPLE_FILES  }                           from '../../modules/local/updater/collect_sample_files'
include { COLLECT_PROJECT_FILES }                           from '../../modules/local/updater/collect_project_files'
include { CREATE_FAIRY_FILE     }                           from '../../modules/local/updater/create_fairy_file'

// ANSI escape code for orange (bright yellow)
def orange = '\033[38;5;208m'
def reset = '\033[0m'

workflow CREATE_INPUT_CHANNELS {
    take:
        indir        // params.indir
        samplesheet  // params.input
        centar       // true when centar is run

    main:
        ch_versions = Channel.empty() // Used to collect the software versions

        //if input directory is passed use it to gather assemblies otherwise use samplesheet
        if (indir != null) {
            indir.view { "CIC_indir: $it" }
//            println("IDGC: ${indir.getClass()}")
            def pattern = params.indir.toString()
            println("Pattern 1 : ${pattern} ${pattern.getClass()}")

            // Create a channel that emits the names of all directories inside the parent directory - we will use this if no samples have a fairy file
            dir_names_ch = indir.map{ parent -> parent.listFiles()  // Load the directory as a Path channel
                    .findAll { it.isDirectory() && it.listFiles().any{ file -> file.name.endsWith('_summaryline.tsv') } }}.flatten().map{dir -> dir.getName().replaceFirst(pattern, '').replaceFirst("/", '')}.collect() // Filter only directories and Filter out excluded directories

            dir_names_ch.view { "Directory names channel: $it" }

            //get list of all samples in the folder - just using the file_integrity file to check that samples that failed are removed from pipeline
            def file_integrity_glob = append_to_path(params.indir.toString(),'*/file_integrity/*_summary.txt')
            // loop through files and identify those that don't have "FAILED" in them and then parse file name and return those ids that pass
            passed_id_ch = Channel.fromPath(file_integrity_glob).collect().map{ it -> get_only_passing_samples(it)}.filter{ it != null }.ifEmpty([])

            passed_id_ch.view { "Passed ID channel: $it" }

            // find the samples that do not have a fairy file and create them
            no_fairy_file_id_ch = passed_id_ch.toList().combine(dir_names_ch.toList()).map{ passed_ids, dir_names -> dir_names - passed_ids}.ifEmpty([])
            no_fairy_file_id_ch.view { "No fairy file ID channel: $it" }
            // print out the samples that did not have a fairy file yet
            no_fairy_file_id_ch.subscribe { result -> if (!result.isEmpty()) { def flat_results = result.flatten().unique().collect() // Check if the channel is empty and print warning for user as this is required for the pipeline
                println("${orange}Warning: There are no files in */file_integrity/*_summary.txt for ${flat_results}. This/These file(s) is required so we will make it.${reset}")}}
            no_fairy_file_id_ch.view { "No fairy file ID channel after subscription: $it" }
            // Fallback to dir_names_ch if passed_id_channel is empty -- need to get all sample ids
            //passed_id_channel = passed_id_ch.concat(dir_names_ch).flatten().unique().collect().toList()
            passed_id_channel = passed_id_ch.concat(dir_names_ch).flatten().flatten().unique().collect().toList()
            passed_id_channel.view { "Passed ID channel after concat: $it" }
            // To make things backwards compatible we need to check if the file_integrity sample is there and if not create it.
            // Collect all ids and combine with indir
            isolates_that_need_file_integrity_ch = indir.map{ it -> get_ids(it) }.flatten().combine(indir).combine(no_fairy_file_id_ch).map{ tuple -> def (old_meta, indir, no_fairy_file_id) = tuple  // Safely destructure the combined tuple
                def meta = [:]
                meta.id = old_meta
                def cleanPath = indir.toString().startsWith('./') ? indir.toString()[2..-1] : indir.toString()
                def cleanerPath = cleanPath.toString().split('/')[-1]
                meta.project_id = cleanerPath
                //def cleaned_path = new File(cleanerPath).getAbsolutePath()
                def cleaned_path = new File(indir.toString()).getAbsolutePath()
                def file_integrity_exists = !(no_fairy_file_id ?: []).contains(meta.id) // Check if the sample is in has no fairy file
                //def file_integrity_exists = no_fairy_file_id.contains(meta.id)
                return [ meta, cleaned_path, file_integrity_exists ]} //.filter{ meta, dir, file_integrity_exists -> file_integrity_exists == false } // Filter out samples that already have a fairy file

            isolates_that_need_file_integrity_ch.view { "Isolates that need file integrity channel: $it" }

            // Now that we have list of samples that need fairy files created make them
            CREATE_FAIRY_FILE (
                isolates_that_need_file_integrity_ch, true
            )
            ch_versions = ch_versions.mix(CREATE_FAIRY_FILE.out.versions)

            // get file_integrity file for MLST updating
            def scaffolds_integrity_glob = append_to_path(params.indir.toString(),'*/file_integrity/*_summary.txt')
           //create file_integrity file channel with meta information -- need to pass to DO_MLST subworkflow
            glob_file_integrity_ch = Channel.fromPath(scaffolds_integrity_glob) // use created regrex to get samples
                    .map{ it -> create_meta(it, "_summary.txt", params.indir.toString(), true)}
                    .combine(passed_id_channel).filter{ meta, file_integrity, passed_id_channel -> passed_id_channel.contains(meta.id)} // check if passed_id_channel contains the id value from the file_integrity file
                    .map{ meta, file_integrity, passed_id_channel -> [meta, file_integrity]} //remove id_channel from output

            //get all integrity files - ones that where already made and ones we just made.
            file_integrity_ch = CREATE_FAIRY_FILE.out.created_fairy_file.collect().ifEmpty([]).combine(glob_file_integrity_ch.ifEmpty([])).flatten().unique().buffer(size:2)

            // loop through files and identify those that don't have "FAILED" in them and then parse file name and return those ids that pass
            all_passed_id_channel = file_integrity_ch.map{meta, dir -> dir}.collect().map{ it -> get_only_passing_samples(it)}.filter { it != null }.toList()

            /////////////////////////// COLLECT ISOLATE LEVEL FILES ///////////////////////////////

            //make relative path full and get 
            directory_ch = all_passed_id_channel.flatten().combine(indir.map{ dir -> 
                def cleanPath = dir.toString().startsWith('./') ? dir.toString()[2..-1] : dir.toString()
                return new File(cleanPath.toString()).getAbsolutePath()})
                .map{ id, dir ->
                        def meta = [:] // create meta array
                        meta.id = id
                        meta.project_id = dir.toString().split('/')[-1]
                        return [meta, dir]} //even though there is only one directory we have to add this so the code works for indir and input

            // Get reads -- tuple -> def (old_meta, indir, passed_ids) = tuple
            def r1_glob = append_to_path(params.indir.toString(),'*/fastp_trimd/*_1.trim.fastq.gz')
            def r2_glob = append_to_path(params.indir.toString(),'*/fastp_trimd/*_2.trim.fastq.gz')
            //create reads channel with meta information
            filtered_r1_reads_ch = Channel.fromPath(r1_glob).map{ it -> create_meta(it, "_1.trim.fastq.gz", params.indir.toString(),false)}
                .combine(all_passed_id_channel).filter{ meta, read_1, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)}
                .map{ meta, read_1, all_passed_id_channel -> [meta, read_1]} //remove all_passed_id_channel from output
            filtered_r2_reads_ch = Channel.fromPath(r2_glob).map{ it -> create_meta(it, "_2.trim.fastq.gz", params.indir.toString(),false)}
                .combine(all_passed_id_channel).filter{ meta, read_2, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)}
                .map{ meta, read_2, all_passed_id_channel -> [meta, read_2]} //remove all_passed_id_channel from output
            // combine reads into one channel
            combined_reads_ch = filtered_r1_reads_ch.join(filtered_r2_reads_ch, by: [0]).map{ meta, read_1, read_2 -> [meta, [read_1, read_2]]}

            // Get scaffolds
            def scaffolds_glob = append_to_path(params.indir.toString(),'*/assembly/*.filtered.scaffolds.fa.gz')
            //create scaffolds channel with meta information -- annoying, but you have to keep this in the brackets instead of having it once outside.
            scaffolds_ch = Channel.fromPath(scaffolds_glob).map{ it -> create_meta(it, ".filtered.scaffolds.fa.gz", params.indir.toString(),false)} // use created regrex to get samples
            // Checking regrex has correct extension
            filtered_scaffolds_ch = scaffolds_ch.map{ it -> check_scaffolds(it) } // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, scaffolds, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} // Filter other channels based on meta.id
                .map{ meta, scaffolds, all_passed_id_channel -> [meta, scaffolds]} //remove all_passed_id_channel from output

            // get .tax files for MLST updating
            def taxa_glob = append_to_path(params.indir.toString(),'*/*.tax')
            //create .tax file channel with meta information 
            filtered_taxonomy_ch = Channel.fromPath(taxa_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, ".tax", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, taxa, all_passed_id_channel-> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, taxa, all_passed_id_channel -> [meta, taxa]} //remove all_passed_id_channel from output

            // get .tax files for MLST updating
            def trimmed_stats_glob = append_to_path(params.indir.toString(),'*/qc_stats/*_trimmed_read_counts.txt')
            //create .tax file channel with meta information 
            filtered_trimmed_stats_ch = Channel.fromPath(trimmed_stats_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_trimmed_read_counts.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, trimmed_stats, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, trimmed_stats, all_passed_id_channel -> [meta, trimmed_stats]} //remove all_passed_id_channel from output

            // get *.top_kraken_hit.txt files for MLST updating
            def kraken_bh_glob = append_to_path(params.indir.toString(),'*/kraken2_trimd/*.kraken2_trimd.top_kraken_hit.txt')
            //create *.top_kraken_hit.txt file channel with meta information 
            filtered_kraken_bh_ch = Channel.fromPath(kraken_bh_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, ".kraken2_trimd.top_kraken_hit.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, kraken_bh, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, kraken_bh, all_passed_id_channel -> [meta, kraken_bh]} //remove all_passed_id_channel from output

            // get .gamma files for gamma HV updating
            def gamma_hv_glob = append_to_path(params.indir.toString(),'*/gamma_hv/*.gamma')
            //create .tax file channel with meta information 
            filtered_gamma_hv_ch = Channel.fromPath(gamma_hv_glob) // use created regrex to get samples
                .map{ it -> create_meta_non_extension(it, params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, gamma_hv, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, gamma_hv, all_passed_id_channel -> [meta, gamma_hv]} //remove all_passed_id_channel from output

            // get .tax files for MLST updating
            def gamma_pf_glob = append_to_path(params.indir.toString(),'*/gamma_pf/*.gamma')
            //create .gamma file channel with meta information 
            filtered_gamma_pf_ch = Channel.fromPath(gamma_pf_glob) // use created regrex to get samples
                .map{ it -> create_meta_non_extension(it, params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, gamma_pf, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, gamma_pf, all_passed_id_channel -> [meta, gamma_pf]} //remove all_passed_id_channel from output

            // get .tax files for MLST updating
            def quast_glob = append_to_path(params.indir.toString(),'*/quast/*_summary.tsv')
            //create .tax file channel with meta information 
            filtered_quast_ch = Channel.fromPath(quast_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_summary.tsv", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, quast_report, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, quast_report, all_passed_id_channel -> [meta, quast_report]} //remove all_passed_id_channel from output

            // get .tax files for MLST updating
            def ani_glob = append_to_path(params.indir.toString(),'*/ANI/*.ani.txt')
            //create .tax file channel with meta information 
            filtered_ani_ch = Channel.fromPath(ani_glob) // use created regrex to get samples
                .map{ it -> create_meta_with_wildcard(it, ".ani.txt", params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, ani, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, ani, all_passed_id_channel -> [meta, ani]} //remove all_passed_id_channel from output

            // get .fastANI.txt
            def fastani_glob = append_to_path(params.indir.toString(),'*/ANI/*.fastANI.txt')
            //create .tax file channel with meta information 
            filtered_ani_best_hit_ch = Channel.fromPath(fastani_glob) // use created regrex to get samples
                .map{ it -> create_meta_with_wildcard(it, ".fastANI.txt", params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, ani_best_hit, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, ani_best_hit, all_passed_id_channel -> [meta, ani_best_hit]} //remove all_passed_id_channel from output

            def assembly_ratio_glob = append_to_path(params.indir.toString(),'*/*_Assembly_ratio_*.txt')
            //create .tax file channel with meta information 
            filtered_assembly_ratio_ch = Channel.fromPath(assembly_ratio_glob) // use created regrex to get samples
                .map{ it -> create_meta_non_extension(it, params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, assembly_ratio, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, assembly_ratio, all_passed_id_channel -> [meta, assembly_ratio]} //remove all_passed_id_channel from output

            // get prokka files
            def prokka_gff_glob = append_to_path(params.indir.toString(),'*/annotation/*.gff')
            //create .gff and .faa files channel with meta information 
            filtered_gff_ch = Channel.fromPath(prokka_gff_glob) // use created regrex to get samples 
                .map{ it -> create_meta_non_extension(it,params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, gff, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, gff, all_passed_id_channel -> [meta, gff]} //remove all_passed_id_channel from output

            def prokka_faa_glob = append_to_path(params.indir.toString(),'*/annotation/*.faa')
            filtered_faa_ch = Channel.fromPath(prokka_faa_glob) // use created regrex to get samples
                .map{ it -> create_meta_non_extension(it, params.indir.toString())}
                .combine(all_passed_id_channel).filter{ meta, faa, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, faa, all_passed_id_channel -> [meta, faa]} //remove all_passed_id_channel from output

            def line_summary_glob = append_to_path(params.indir.toString(),'*/*_summaryline.tsv')
            line_summary_ch = Channel.fromPath(line_summary_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_summaryline.tsv", params.indir.toString(),false)} // create meta for sample

            //////////////////////////// FOR CENTAR ///////////////////////////////////////

            // get files for MLST updating 
            def combined_mlst_glob = append_to_path(params.indir.toString(),'*/mlst/*_combined.tsv')

            //create .tax file channel with meta information 
            filtered_combined_mlst_ch = Channel.fromPath(combined_mlst_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_combined.tsv", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, combined_mlst, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples - keep those in all_passed_id_channel
                .map{ meta, combined_mlst, all_passed_id_channel -> [meta, combined_mlst]} //remove all_passed_id_channel from output

            /////////////////////////// COLLECT SPECIES SPECIFIC FILES /////////////////////////////

            // get files for MLST updating 
            def shigapass_glob = append_to_path(params.indir.toString(),'*/ANI/*_ShigaPass_summary.csv')

            //collect .tax file channel with meta information 
            shigapass_files_ch = Channel.fromPath(shigapass_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_ShigaPass_summary.csv", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, shiapass_files, all_passed_id -> all_passed_id.contains(meta.id)} //filtering out failured samples - keep those in all_passed_id_channel
                .map{ meta, shiapass_files, all_passed_id -> [meta, shiapass_files]} //remove all_passed_id_channel from output'

            // get CENTAR files
            def centar_glob = append_to_path(params.indir.toString(),'*/CENTAR/*_centar_output.tsv')

            //collect .tax file channel with meta information 
            centar_files_ch = Channel.fromPath(centar_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_centar_output.tsv", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, centar_files, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples - keep those in all_passed_id_channel
                .map{ meta, centar_files, all_passed_id_channel -> [meta, centar_files]}.ifEmpty( [[id: "", project_id: ""], []] ) //remove all_passed_id_channel from output

            /////////////////////////// COLLECT README FOR UPDATER ///////////////////////////////
            def readme_glob = append_to_path(params.indir.toString(),'*/*_updater_log.tsv')

            /*Channel.fromPath(readme_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_updater_log.tsv", params.indir.toString(), false)}.ifEmpty( [[id: "", project_id: ""], []] ) // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, readme_files, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)}.ifEmpty( [[id: "", project_id: ""], [], []] )
                .map{ meta, readme_files, all_passed_id_channel -> [meta, readme_files]}.view()*/

            readme_files_ch = Channel.fromPath(readme_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_updater_log.tsv", params.indir.toString(), false)}.ifEmpty( [[id: "", project_id: ""], []] )  // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, readme_files, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)}.ifEmpty( [[id: "", project_id: ""], [], []] ) //filtering out failed samples - keep those in all_passed_id_channel
                .map{ meta, readme_files, all_passed_id_channel -> [meta, readme_files]} //remove all_passed_id_channel from output

            /////////////////////////// COLLECT PROJECT LEVEL FILES ///////////////////////////////

            def synopsis_glob = append_to_path(params.indir.toString(),'*/*.synopsis')
            synopsis_ch = Channel.fromPath(synopsis_glob) // use created regrex to get samples
                .map{ it -> create_meta_non_extension(it, params.indir.toString())} // create meta for sample and adding group ID
            //filtering out failured samples
            filtered_synopsis_ch = synopsis_ch.filter{ meta, synopsis -> all_passed_id_channel.contains(meta.id)}

            // use created regrex to get samples
            def griphin_excel_glob = append_to_path(params.indir.toString(),'*_GRiPHin_Summary.xlsx')
            all_griphin_excel_ch = all_passed_id_channel.flatten().combine(Channel.fromPath(griphin_excel_glob)).map{ it -> create_groups_and_id(it, params.indir.toString())}
                //.map{ it -> modifiedFileChannel(it, "_GRiPHin_Summary","_old_GRiPHin") }
            //create regrex, get files in dir, add in meta information, change name
            def griphin_tsv_glob = append_to_path(params.indir.toString(),'*_GRiPHin_Summary.tsv')
            all_griphin_tsv_ch = all_passed_id_channel.flatten().combine(Channel.fromPath(griphin_tsv_glob)).map{ it -> create_groups_and_id(it, params.indir.toString())} 
                //.map{ it -> modifiedFileChannel(it, "_GRiPHin_Summary","_old_GRiPHin") }
            def phoenix_tsv_glob = append_to_path(params.indir.toString(),'Phoenix_Summary.tsv')
            all_phoenix_tsv_ch = all_passed_id_channel.flatten().combine(Channel.fromPath(phoenix_tsv_glob)).map{ it -> create_groups_id_and_busco(it, params.indir.toString())}
                //.map{ it -> modifiedFileChannel(it, "_Summary","_Summary_Old") }
            def pipeline_info_glob = append_to_path(params.indir.toString(),'pipeline_info/software_versions.yml')
            all_pipeline_info_ch = all_passed_id_channel.flatten().combine(Channel.fromPath(pipeline_info_glob)).map{ it -> create_groups_and_id(it, params.indir.toString())} // use created regrex to get samples

            // this if/else is only here to make sure the output goes to the correct output folder as its different for each in the modules.config.
            if (centar == true) {
                //get valid samplesheet for griphin step in cdc_scaffolds
                CENTAR_CREATE_SAMPLESHEET (
                    indir
                )
                ch_versions = ch_versions.mix(CENTAR_CREATE_SAMPLESHEET.out.versions)
                
                valid_samplesheet = CENTAR_CREATE_SAMPLESHEET.out.samplesheet
            } else {
                //get valid samplesheet for griphin step in cdc_scaffolds
                CREATE_SAMPLESHEET (
                    indir
                )
                ch_versions = ch_versions.mix(CREATE_SAMPLESHEET.out.versions)

                valid_samplesheet = CREATE_SAMPLESHEET.out.samplesheet
            }

            // combining all summary files into one channel
            summary_files_ch = all_griphin_excel_ch.join(all_griphin_tsv_ch.map{meta, griphin_tsv   -> [[id:meta.id, project_id:meta.project_id], griphin_tsv]},   by: [[0][0],[0][1]])
                                    .join(all_phoenix_tsv_ch.map{               meta, phoenix_tsv   -> [[id:meta.id, project_id:meta.project_id], phoenix_tsv]},   by: [[0][0],[0][1]])
                                    .join(all_pipeline_info_ch.map{             meta, pipeline_info -> [[id:meta.id, project_id:meta.project_id], pipeline_info]}, by: [[0][0],[0][1]])

            // pulling all the necessary project level files into channels - need to do this for the name change.
            COLLECT_PROJECT_FILES (
                summary_files_ch, false
            )
            ch_versions = ch_versions.mix(COLLECT_PROJECT_FILES.out.versions)

            griphin_excel_ch = COLLECT_PROJECT_FILES.out.griphin_excel
            griphin_tsv_ch = COLLECT_PROJECT_FILES.out.griphin_tsv
            phoenix_tsv_ch = COLLECT_PROJECT_FILES.out.phoenix_tsv.map{it -> add_entry_meta(it)}
            pipeline_info_ch = COLLECT_PROJECT_FILES.out.software_versions_file

        } else if (samplesheet != null) {

            if (centar == true) { // this if/else is only here to make sure the output goes to the correct output folder as its different for each in the modules.config.
                // if a samplesheet was passed then use that to create the channel
                CENTAR_SAMPLESHEET_CHECK (
                    samplesheet, false, false, true
                )
                ch_versions = ch_versions.mix(CENTAR_SAMPLESHEET_CHECK.out.versions)

                samplesheet = CENTAR_SAMPLESHEET_CHECK.out.csv
            } else {
                // if a samplesheet was passed then use that to create the channel
                SAMPLESHEET_CHECK (
                    samplesheet, false, false, true
                )
                ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

                samplesheet = SAMPLESHEET_CHECK.out.csv
            }

            // To make things backwards compatible we need to check if the file_integrity sample is there and if not create it.
            file_integrity_exists = samplesheet.splitCsv( header:true, sep:',' ).map{ it -> check_file_integrity(it) }

            CREATE_FAIRY_FILE (
                file_integrity_exists, false
            )
            ch_versions = ch_versions.mix(CREATE_FAIRY_FILE.out.versions)

            directory_ch = samplesheet.splitCsv( header:true, sep:',' ).map{ it -> create_dir_channels(it) }

            //adding meta.id to end of dir - otherwise too many files are copied and it takes forever. 
            sample_directory_ch = samplesheet.splitCsv( header:true, sep:',' ).map{ it -> create_sample_dir_channels(it) }

            // pulling all the necessary sample level files into channels
            COLLECT_SAMPLE_FILES (
                sample_directory_ch
            )
            ch_versions = ch_versions.mix(COLLECT_SAMPLE_FILES.out.versions)

            //collect all fairy files and then recreate meta groups with flatten and buffer.
            file_integrity_ch = CREATE_FAIRY_FILE.out.created_fairy_file.collect().ifEmpty([]).combine(COLLECT_SAMPLE_FILES.out.fairy_summary.collect().ifEmpty([])).flatten().buffer(size:2)

            //combine reads to get into one channel
            combined_reads_ch = COLLECT_SAMPLE_FILES.out.read1.join(COLLECT_SAMPLE_FILES.out.read2, by: [0]).map{ meta, read1, read2 -> [meta, [read1, read2]]}
            // get other files
            filtered_scaffolds_ch = COLLECT_SAMPLE_FILES.out.scaffolds
            filtered_gff_ch = COLLECT_SAMPLE_FILES.out.gff
            filtered_faa_ch = COLLECT_SAMPLE_FILES.out.faa
            line_summary_ch = COLLECT_SAMPLE_FILES.out.summary_line
            filtered_synopsis_ch = COLLECT_SAMPLE_FILES.out.synopsis
            filtered_taxonomy_ch = COLLECT_SAMPLE_FILES.out.tax
            filtered_gamma_pf_ch = COLLECT_SAMPLE_FILES.out.gamma_pf
            filtered_gamma_hv_ch = COLLECT_SAMPLE_FILES.out.gamma_hv
            filtered_assembly_ratio_ch = COLLECT_SAMPLE_FILES.out.assembly_ratio
            filtered_kraken_bh_ch = COLLECT_SAMPLE_FILES.out.kraken_bh
            filtered_trimmed_stats_ch = COLLECT_SAMPLE_FILES.out.trimmed_stats
            filtered_quast_ch = COLLECT_SAMPLE_FILES.out.quast_report
            filtered_ani_ch = COLLECT_SAMPLE_FILES.out.ani
            filtered_ani_best_hit_ch = COLLECT_SAMPLE_FILES.out.ani_best_hit
            filtered_combined_mlst_ch = COLLECT_SAMPLE_FILES.out.combined_mlst

            //species specific files
            shigapass_files_ch = COLLECT_SAMPLE_FILES.out.shigapass_output
            centar_files_ch = COLLECT_SAMPLE_FILES.out.centar_output
            //readme files
            readme_files_ch = COLLECT_SAMPLE_FILES.out.readme

            summary_files_ch = samplesheet.splitCsv( header:true, sep:',' ).map{ it -> create_summary_files_channels(it) }

            // pulling all the necessary project level files into channels
            COLLECT_PROJECT_FILES (
                summary_files_ch, true
            )
            ch_versions = ch_versions.mix(COLLECT_PROJECT_FILES.out.versions)

            griphin_excel_ch = COLLECT_PROJECT_FILES.out.griphin_excel
            griphin_tsv_ch = COLLECT_PROJECT_FILES.out.griphin_tsv
            phoenix_tsv_ch = COLLECT_PROJECT_FILES.out.phoenix_tsv.map{it -> add_entry_meta(it)}
            pipeline_info_ch = COLLECT_PROJECT_FILES.out.software_versions_file

            valid_samplesheet = samplesheet // still need to check this file

        } else {
            exit 1, 'You need EITHER an input directory samplesheet (using --input) or a single sample directory (using --indir)!' 
        }

    emit:
        //project level summary files
        griphin_excel_ch   = griphin_excel_ch
        griphin_tsv_ch     = griphin_tsv_ch
        phoenix_tsv_ch     = phoenix_tsv_ch
        pipeline_info_ch   = pipeline_info_ch
        directory_ch       = directory_ch
        valid_samplesheet  = valid_samplesheet
        versions           = ch_versions
        pipeline_info      = pipeline_info_ch

        //species specific files
        centar             = centar_files_ch
        shigapass          = shigapass_files_ch
        //updater
        readme             = readme_files_ch

        // sample specific files
        filtered_scaffolds = filtered_scaffolds_ch      // channel: [ meta, [ scaffolds_file ] ]
        reads              = combined_reads_ch
        taxonomy           = filtered_taxonomy_ch
        prokka_gff         = filtered_gff_ch
        prokka_faa         = filtered_faa_ch
        fairy_outcome      = file_integrity_ch
        line_summary       = line_summary_ch // need non-filtered to make summary files will all samples in project folder
        synopsis           = filtered_synopsis_ch
        ani                = filtered_ani_ch
        ani_best_hit       = filtered_ani_best_hit_ch
        gamma_pf           = filtered_gamma_pf_ch
        gamma_hv           = filtered_gamma_hv_ch
        assembly_ratio     = filtered_assembly_ratio_ch
        k2_bh_summary      = filtered_kraken_bh_ch
        fastp_total_qc     = filtered_trimmed_stats_ch
        quast_report       = filtered_quast_ch
        combined_mlst      = filtered_combined_mlst_ch // for centar entry

}

/*========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/

def check_file_integrity(LinkedHashMap row) {
    def meta = [:] // create meta array
    meta.id = row.sample
    //meta.project_id = row.directory.toString().split('/')[-2]
    //meta.project_id = row.directory
    def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
    meta.project_id = new File(clean_path).getParent()
    def pattern = "*_summary.txt"
    // Convert the wildcard pattern to regex: "*_summary.txt" to ".*_summary\.txt"
    def regexPattern = pattern.replace("*", ".*").replace("?", ".")
    File dir = new File(clean_path + "/file_integrity/")
    // List files matching the regex pattern
    def files = dir.listFiles { file -> file.name ==~ /${regexPattern}/ }
    if (files && files.length > 0) {
        files.each { file ->
            def lines = file.readLines()
            if (lines.size() != 5) {
                exit 1, "ERROR: File '${file.name}' in '${file.parent}' has ${lines.size()} lines instead of 5, this will cause errors downstream, please fix and rerun."
            }
        }
        return [ meta, clean_path, true ]
    } else {
        return [ meta, clean_path, false ]
    }
}

def get_ids(dir) {
    // Initialize an empty list to store directory names
    List<String> dirNames = []

    // Check if the path is indeed a directory
    if (dir.isDirectory()) {
        // Collect only first-level subdirectories
        dir.eachFile { file ->
            if (file.isDirectory()) {
                // Add directory name if it's not in the exclude list
                if (file.listFiles().any { it.name.endsWith('_summaryline.tsv')}) {
                    dirNames << file.name
                }
            }
        }
    return dirNames
    }
}

def modifiedFileChannel(input_ch, old_string, new_string) {
    def newFileName = input_ch[1].getName().replaceAll(old_string, new_string)
    def newFilePath = input_ch[1].getParent() ? input_ch[1].getParent().resolve(newFileName) : newFileName
    return [ input_ch[0], newFilePath ]
}

// Function to get list of [ meta, [ directory ] ]
def create_summary_files_channels(LinkedHashMap row) {
    def meta = [:] // create meta array
    meta.id = row.sample
    def project_id = row.directory.toString().split('/')[-2]
    //meta.project_id = row.directory
    def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
    def cleaned_path = new File(clean_path).getParent()
    meta.project_id = cleaned_path
    def software_versions = cleaned_path + "/pipeline_info/"
    def griphin_summary_tsv = cleaned_path + "/" + project_id + "_GRiPHin_Summary.tsv"
    def griphin_summary_excel = cleaned_path + "/" + project_id +  "_GRiPHin_Summary.xlsx"
    def phx_summary = cleaned_path + "/Phoenix_Summary.tsv"
    return [ meta, griphin_summary_excel, griphin_summary_tsv, phx_summary, software_versions ]
}

// Function to get list of [ meta, [ directory ] ]
def create_dir_channels(LinkedHashMap row) {
    def meta = [:] // create meta array
    meta.id = row.sample
    //meta.project_id = row.directory.toString().split('/')[-2]
    //meta.project_id = row.directory
    def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
    def cleaned_path = new File(clean_path).getParent()
    meta.project_id = cleaned_path
    return [ meta, cleaned_path ]
}

def add_entry_meta(input_ch){
    """samples needed to be grouped by their project directory for editing summary files."""
    def meta = [:] // create meta array
    meta.id = input_ch[0].id
    meta.project_id = input_ch[0].project_id
    // Use file object to read the first line of the file
    def file = input_ch[1]
    def firstLine = file.text.split('\n')[0]
    meta.entry = firstLine.contains('BUSCO')
    return [meta, input_ch[1]]
}

def create_groups_id_and_busco(input_ch, project_folder){
    """samples needed to be grouped by their project directory for editing summary files."""
    def meta = [:] // create meta array
    meta.id = input_ch[0]
    meta.project_id = project_folder.toString().split('/')[-1]
    def firstLine = input_ch[1].text.split('\n')[0]
    meta.entry = firstLine.contains('BUSCO')
    return [meta, input_ch[1]]
}

// Function to get list of [ meta, [ directory ] ]
def create_sample_dir_channels(LinkedHashMap row) {
    def meta = [:] // create meta array
    meta.id = row.sample
    //meta.project_id = row.directory.toString().split('/')[-2]
    //meta.project_id = row.directory
    def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
    meta.project_id = new File(clean_path).getParent()
    return [ meta, clean_path ]
}

def create_groups_and_id(input_ch, project_folder){
    """samples needed to be grouped by their project directory for editing summary files."""
    def meta = [:] // create meta array
    meta.id = input_ch[0]
    meta.project_id = project_folder.toString().split('/')[-1]
    return [meta, input_ch[1]]
}

def create_groups(input_ch, project_folder){
    """samples needed to be grouped by their project directory for editing summary files."""
    def meta = [:] // create meta array
    meta.id = input_ch[0]
    meta.project_id = project_folder.toString()
    return [meta, project_folder]
}

// Function to filter files and parse their names
def get_only_passing_samples(filePaths) {
    def passingFiles = [] // Initialize an empty list to store passing file names
    filePaths.each { filePath ->
        // Exclude file paths containing ".nextflow"
        if (filePath.toString().contains('.nextflow')) {
            println("Skipping: ${filePath} (contains .nextflow dir).")
            return  // Skip this iteration
        }
        def fileObj = new File(filePath.toString())
        if (!fileObj.text.contains('FAILED')) {
            // Get the file name by removing the path and extension
            def fileName = filePath.getFileName().toString().split('/')[-1]
            // Remove everything after the regex '_*_summary.txt'
            def passing_file1 = fileName.replaceAll(/_scaffolds_summary\.txt/, '').replaceAll(/_rawstats_summary\.txt/, '').replaceAll(/_corruption_summary\.txt/, '').replaceAll(/_trimstats_summary\.txt/, '')
            def passing_file2 = passing_file1.replaceAll(/_summary\.txt/, '')
            passingFiles << passing_file2 // Add the passing file name to the list
        }
    }
    return passingFiles
}

def append_to_path(full_path, string) {
    if (full_path.toString().endsWith('/')) {
        return full_path.toString() + string
    } else {
        return full_path.toString() + '/' + string
    }
}

def create_meta_with_wildcard(sample, file_extension, indir){
    '''Creating meta: [[id:sample1], $PATH/sample1_REFSEQ_20240124.fastANI.txt]'''
    def meta = [:] // create meta array
    meta.id = sample.getName().replaceAll(file_extension, "").split('_REFSEQ_\\d{8}')[0] // get file name without extention
    meta.project_id = indir.toString().split('/')[-1]
    return [ meta, sample ]
}

def create_meta(sample, file_extension, indir, extra_check){
    '''Creating meta: [[id:sample1], $PATH/sample1.filtered.scaffolds.fa.gz]'''
    def meta = [:] // create meta array
    meta.id = sample.getName().replaceAll(file_extension, "") // get file name without extention
    if (extra_check==true){
        full_ext1 =  "_scaffolds" + file_extension
        full_ext2 = "_corruption" + file_extension
        full_ext3 = "_trimstats" + file_extension
        full_ext4 = "_rawstats" + file_extension
        meta.id = sample.getName().replaceAll(full_ext1, "").replaceAll(full_ext2, "").replaceAll(full_ext3, "").replaceAll(full_ext4, "").replaceAll(file_extension, "")
    } else {
        meta.id = sample.getName().replaceAll(file_extension, "") // get file name without extention
    }
    meta.project_id = indir.toString().split('/')[-1]
    return [ meta, sample ]
}

def create_meta_non_extension(sample, indir){
    '''Creating meta: [[id:sample1], $PATH/sample1.filtered.scaffolds.fa.gz]'''
    def meta = [:] // create meta array
    meta.id = sample.getSimpleName()//.split('_')[0] // get the last string after the last backslash
    // Regex patterns to match _HyperVirulence_YYYYMMDD and _PF-Replicons_YYYYMMDD
    def hyperVirulencePattern = /_HyperVirulence_\d{8}/
    def pfRepliconsPattern = /_PF-Replicons_\d{8}/
    def assemblyratioPattern = /_Assembly_ratio_\d{8}/

    // Check if the id contains either of the patterns
    if (meta.id =~ hyperVirulencePattern || meta.id =~ pfRepliconsPattern || meta.id =~ assemblyratioPattern) {
        // Remove the pattern if it matches
        meta.id = meta.id.replaceAll(hyperVirulencePattern, '').replaceAll(pfRepliconsPattern, '').replaceAll(assemblyratioPattern, '').trim()} // Trim any trailing or leading spaces

    meta.project_id = indir.toString().split('/')[-1]
    return [ meta, sample ]
}

def check_scaffolds(scaffold_channel) {
    if (scaffold_channel[1].toString().endsWith(".fasta.gz") || scaffold_channel[1].toString().endsWith(".fa.gz") ) {
        return scaffold_channel
    } else {
        exit 1, "ERROR: No scaffolds found in '${scaffold_channel[1].toString()}'. Either your scaffold regrex is off (scaffolds files should end in either '.fa.gz' or '.fasta.gz') or the directory provided doesn't contain assembly files." 
    }
}
