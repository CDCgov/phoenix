//
// workflow handles taking in either a samplesheet or directory and creates correct channels for scaffolds entry point
//
// for centar entry
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

// Add this outside your workflow block
def check_update(meta, file, db, type, mode, needsUpdate = false) {
    return InputChannelUtils.previous_updater_check(meta, file, db, type, mode, needsUpdate)
}

workflow CREATE_INPUT_CHANNELS {
    take:
        indir        // params.indir
        samplesheet  // params.input
        centar       // true when centar is run

    main:

        ch_versions = Channel.empty() // Used to collect the software versions

        //if input directory is passed use it to gather assemblies otherwise use samplesheet
        if (indir != null) {
            def pattern = params.indir.toString()

            // Create a channel that emits the names of all directories inside the parent directory - we will use this if no samples have a fairy file
            dir_names_ch = indir.map{ parent -> parent.listFiles()  // Load the directory as a Path channel
                    .findAll { it.isDirectory() && it.listFiles().any{ file -> file.name.endsWith('_summaryline.tsv') } }}.flatten().map{dir -> dir.getName().replaceFirst(pattern, '').replaceFirst("/", '')}.collect() // Filter only directories and Filter out excluded directories

            //get list of all samples in the folder - just using the file_integrity file to check that samples that failed are removed from pipeline
            def file_integrity_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/file_integrity/*_summary.txt')
            // create a channel with the ids of the samples that failed the file integrity checks to print a warning for the user
            failed_fairy_ids_ch = Channel.fromPath(file_integrity_glob).collect().map{ it -> InputChannelUtils.get_failed_samples(it)}.filter{ it != null }.ifEmpty([])
            failed_fairy_ids_ch.subscribe { result -> if (!result.isEmpty()) { def flat_results = result.flatten().unique().collect() // Check if the channel is empty and print warning for user as this is required for the pipeline
                println("${orange}Warning: The following files failed file integrity checks by phoenix and will not be included in this analysis ${flat_results}. ${reset}")}}

            // loop through files and identify those that don't have "FAILED" in them and then parse file name and return those ids that pass
            passed_id_ch = Channel.fromPath(file_integrity_glob).collect().map{ it -> InputChannelUtils.get_only_passing_samples(it)}.filter{ it != null }.ifEmpty([])

            // find the samples that do not have a fairy file and create them
            no_fairy_file_id_ch = passed_id_ch.toList().combine(dir_names_ch.toList()).combine(failed_fairy_ids_ch.toList()).map{ passed_ids, dir_names, failed_ids -> dir_names - failed_ids - passed_ids}.ifEmpty([])

            // print out the samples that did not have a fairy file yet - we will also add in files that have the v2.1.1 style of fairy file
            old_style_fairy_file_ch = Channel.fromPath(file_integrity_glob).collect().map{ it -> InputChannelUtils.get_old_fairy_samples(it)}.filter{ it != null }.ifEmpty([])
            no_fairy_file_id_ch.combine(old_style_fairy_file_ch.toList()).subscribe { result -> if (!result.isEmpty()) { def flat_results = result.flatten().unique().collect() // Check if the channel is empty and print warning for user as this is required for the pipeline
                println("${orange}Warning: There are no files in */file_integrity/*_summary.txt for ${flat_results}. This/These file(s) is required so we will make it.${reset}")}}

            // Fallback to dir_names_ch if passed_id_channel is empty -- need to get all sample ids
            passed_id_channel = passed_id_ch.concat(dir_names_ch).flatten().flatten().unique().collect().toList()

            // To make things backwards compatible we need to check if the file_integrity sample is there and if not create it.
            // Collect all ids and combine with indir -- issue here!!
            isolates_that_need_file_integrity_ch = indir.map{ it -> InputChannelUtils.get_ids(it) }.flatten().combine(indir).combine(no_fairy_file_id_ch.toList()).map{ tuple -> def (old_meta, indir, no_fairy_file_id) = tuple  // Safely destructure the combined tuple
                def meta = [:]
                meta.id = old_meta
                def cleanPath = indir.toString().startsWith('./') ? indir.toString()[2..-1] : indir.toString()
                def cleanerPath = cleanPath.toString().split('/')[-1]
                meta.project_id = cleanerPath
                def cleaned_path = new File(indir.toString()).getAbsolutePath()
                def file_integrity_exists = !(no_fairy_file_id ?: []).contains(meta.id) // Check if the sample is in has no fairy file
                return [ meta, cleaned_path, file_integrity_exists ]}.filter{ meta, dir, file_integrity_exists -> file_integrity_exists == false } // Filter out samples that already have a fairy file

            // Now that we have list of samples that need fairy files created make them
            CREATE_FAIRY_FILE (
                isolates_that_need_file_integrity_ch, true
            )
            ch_versions = ch_versions.mix(CREATE_FAIRY_FILE.out.versions)

            // get file_integrity file for MLST updating
            def scaffolds_integrity_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/file_integrity/*_summary.txt')
           //create file_integrity file channel with meta information -- need to pass to DO_MLST subworkflow
            glob_file_integrity_ch = Channel.fromPath(scaffolds_integrity_glob) // use created regrex to get samples
                    .map{ it -> InputChannelUtils.create_meta(it, "_summary.txt", params.indir.toString(), true)}
                    .combine(passed_id_channel).filter{ meta, file_integrity, passed_id_channel -> passed_id_channel.contains(meta.id)} // check if passed_id_channel contains the id value from the file_integrity file
                    .map{ meta, file_integrity, passed_id_channel -> [meta, file_integrity]} //remove id_channel from output

            //get all integrity files - ones that where already made and ones we just made.
            file_integrity_ch = CREATE_FAIRY_FILE.out.created_fairy_file.collect().ifEmpty([]).combine(glob_file_integrity_ch.ifEmpty([])).flatten().unique().buffer(size:2)

            // loop through files and identify those that don't have "FAILED" in them and then parse file name and return those ids that pass
            all_passed_id_channel = file_integrity_ch.map{meta, dir -> dir}.collect().map{ it -> InputChannelUtils.get_only_passing_samples(it)}.filter { it != null }.toList()

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

            // Get reads
            def r1_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/fastp_trimd/*_1.trim.fastq.gz')
            def r2_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/fastp_trimd/*_2.trim.fastq.gz')
            //create reads channel with meta information
            filtered_r1_reads_ch = Channel.fromPath(r1_glob).map{ it -> InputChannelUtils.create_meta(it, "_1.trim.fastq.gz", params.indir.toString(),false)}
                .combine(all_passed_id_channel).filter{ meta, read_1, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)}
                .map{ meta, read_1, all_passed_id_channel -> [meta, read_1]} //remove all_passed_id_channel from output
            filtered_r2_reads_ch = Channel.fromPath(r2_glob).map{ it -> InputChannelUtils.create_meta(it, "_2.trim.fastq.gz", params.indir.toString(),false)}
                .combine(all_passed_id_channel).filter{ meta, read_2, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)}
                .map{ meta, read_2, all_passed_id_channel -> [meta, read_2]} //remove all_passed_id_channel from output
            // combine reads into one channel
            combined_reads_ch = filtered_r1_reads_ch.join(filtered_r2_reads_ch, by: [0]).map{ meta, read_1, read_2 -> [meta, [read_1, read_2]]}

            // Get scaffolds
            def scaffolds_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/assembly/*.filtered.scaffolds.fa.gz')
            //create scaffolds channel with meta information
            scaffolds_ch = Channel.fromPath(scaffolds_glob).map{ it -> InputChannelUtils.create_meta(it, ".filtered.scaffolds.fa.gz", params.indir.toString(),false)} // use created regrex to get samples
            // Checking regrex has correct extension
            filtered_scaffolds_ch = scaffolds_ch.map{ it -> InputChannelUtils.check_scaffolds(it) } // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, scaffolds, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} // Filter other channels based on meta.id
                .map{ meta, scaffolds, all_passed_id_channel -> [meta, scaffolds]} //remove all_passed_id_channel from output

            // Get renamed scaffolds
            def renamed_scaffolds_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/assembly/*.renamed.scaffolds.fa.gz')
            //create scaffolds channel with meta information
            renamed_scaffolds_ch = Channel.fromPath(renamed_scaffolds_glob).map{ it -> InputChannelUtils.create_meta(it, ".renamed.scaffolds.fa.gz", params.indir.toString(),false)} // use created regrex to get samples
            // Checking regrex has correct extension
            filtered_renamed_scaffolds_ch = renamed_scaffolds_ch.map{ it -> InputChannelUtils.check_scaffolds(it) } // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, renamed_scaffolds, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} // Filter other channels based on meta.id
                .map{ meta, renamed_scaffolds, all_passed_id_channel -> [meta, renamed_scaffolds]} //remove all_passed_id_channel from output

            // get .tax files for MLST updating
            def taxa_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/*.tax')
            //create .tax file channel with meta information
            filtered_taxonomy_ch = Channel.fromPath(taxa_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".tax", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, taxa, all_passed_id_channel-> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, taxa, all_passed_id_channel -> [meta, taxa]} //remove all_passed_id_channel from output

            // get _raw_read_counts.txt files
            def raw_stats_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/raw_stats/*_raw_read_counts.txt')
            //create .tax file channel with meta information
            filtered_raw_stats_ch = Channel.fromPath(raw_stats_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_raw_read_counts.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, raw_stats, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, raw_stats, all_passed_id_channel -> [meta, raw_stats]} //remove all_passed_id_channel from output

            // get .tax files for MLST updating
            def trimmed_stats_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/qc_stats/*_trimmed_read_counts.txt')
            //create .tax file channel with meta information
            filtered_trimmed_stats_ch = Channel.fromPath(trimmed_stats_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_trimmed_read_counts.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, trimmed_stats, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, trimmed_stats, all_passed_id_channel -> [meta, trimmed_stats]} //remove all_passed_id_channel from output

            // get *.top_kraken_hit.txt files for MLST updating
            def trimd_kraken_bh_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_trimd/*.kraken2_trimd.top_kraken_hit.txt')
            //create *.top_kraken_hit.txt file channel with meta information
            filtered_trimd_kraken_bh_ch = Channel.fromPath(trimd_kraken_bh_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".kraken2_trimd.top_kraken_hit.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, kraken_bh, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, kraken_bh, all_passed_id_channel -> [meta, kraken_bh]} //remove all_passed_id_channel from output

            // get *.summary.txt files
            def trimd_kraken_report_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_trimd/*.kraken2_trimd.summary.txt')
            //create *.summary.txt file channel with meta information
            filtered_trimd_kraken_report_ch = Channel.fromPath(trimd_kraken_report_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".kraken2_trimd.summary.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, kraken_report, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, kraken_report, all_passed_id_channel -> [meta, kraken_report]} //remove all_passed_id_channel from output

            // get *_trimd.html files
            def trimd_kraken_krona_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_trimd/krona/*_trimd.html')
            //create *_trimd.html file channel with meta information
            filtered_trimd_krona_ch = Channel.fromPath(trimd_kraken_krona_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_trimd.html", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, krona, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, krona, all_passed_id_channel -> [meta, krona]} //remove all_passed_id_channel from output

            // get *.top_kraken_hit.txt
            def asmbld_kraken_bh_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_asmbld/*.kraken2_asmbld.top_kraken_hit.txt')
            //create *.top_kraken_hit.txt file channel with meta information
            filtered_asmbld_kraken_bh_ch = Channel.fromPath(asmbld_kraken_bh_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".kraken2_asmbld.top_kraken_hit.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, kraken_bh, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, kraken_bh, all_passed_id_channel -> [meta, kraken_bh]} //remove all_passed_id_channel from output

            // get *.summary.txt files
            def asmbld_kraken_report_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_asmbld/*.kraken2_asmbld.summary.txt')
            //create *.summary.txt file channel with meta information
            filtered_asmbld_kraken_report_ch = Channel.fromPath(asmbld_kraken_report_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".kraken2_asmbld.summary.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, kraken_report, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, kraken_report, all_passed_id_channel -> [meta, kraken_report]} //remove all_passed_id_channel from output

            // get *_asmbld.html files
            def asmbld_kraken_krona_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_asmbld/krona/*_asmbld.html')
            //create *_asmbld.html file channel with meta information
            filtered_asmbld_krona_ch = Channel.fromPath(asmbld_kraken_krona_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_asmbld.html", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, krona, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, krona, all_passed_id_channel -> [meta, krona]} //remove all_passed_id_channel from output

            // get short_summary.specific.*.<sample_id>.filtered.scaffolds.fa.txt
            def busco_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/BUSCO/short_summary.specific.*.filtered.scaffolds.fa.txt')
            Channel.fromPath(busco_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".filtered.scaffolds.fa.txt", params.indir.toString(),false)}
            //create short_summary.specific.*.<sample_id>.filtered.scaffolds.fa.txt file channel with meta information
            filtered_busco_short_summary_ch = Channel.fromPath(busco_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".filtered.scaffolds.fa.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, short_summary, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, short_summary, all_passed_id_channel -> [meta, short_summary]} //remove all_passed_id_channel from output

            // get *.top_kraken_hit.txt
            def wtasmbld_kraken_bh_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_asmbld_weighted/*.kraken2_wtasmbld.top_kraken_hit.txt')
            //create *.top_kraken_hit.txt file channel with meta information
            filtered_wtasmbld_kraken_bh_ch = Channel.fromPath(wtasmbld_kraken_bh_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".kraken2_wtasmbld.top_kraken_hit.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, kraken_bh, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, kraken_bh, all_passed_id_channel -> [meta, kraken_bh]} //remove all_passed_id_channel from output

            // get *.summary.txt files
            def wtasmbld_kraken_report_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_asmbld_weighted/*.kraken2_wtasmbld.summary.txt')
            //create *.summary.txt file channel with meta information
            filtered_wtasmbld_kraken_report_ch = Channel.fromPath(wtasmbld_kraken_report_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, ".kraken2_wtasmbld.summary.txt", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, kraken_report, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, kraken_report, all_passed_id_channel -> [meta, kraken_report]} //remove all_passed_id_channel from output

            // get *_wtasmbld.html files
            def wtasmbld_kraken_krona_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/kraken2_asmbld_weighted/krona/*_wtasmbld.html')
            //create *_wtasmbld.html file channel with meta information
            filtered_wtasmbld_krona_ch = Channel.fromPath(wtasmbld_kraken_krona_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_wtasmbld.html", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, krona, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, krona, all_passed_id_channel -> [meta, krona]} //remove all_passed_id_channel from output

            // get quast files for MLST updating
            def quast_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/quast/*_summary.tsv')
            //create .tax file channel with meta information
            filtered_quast_ch = Channel.fromPath(quast_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_summary.tsv", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, quast_report, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, quast_report, all_passed_id_channel -> [meta, quast_report]} //remove all_passed_id_channel from output

            // get .tax files for MLST updating
            def ani_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/ANI/*.ani.txt')
            //create .tax file channel with meta information
            filtered_ani_ch = Channel.fromPath(ani_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta_with_wildcard(it, ".ani.txt", params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, ani, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, ani, all_passed_id_channel -> [meta, ani]} //remove all_passed_id_channel from output

            // get .fastANI.txt
            def fastani_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/ANI/*.fastANI.txt')
            filtered_ani_best_hit_ch = Channel.fromPath(fastani_glob)
                .map{ it -> InputChannelUtils.create_meta_with_wildcard(it, ".fastANI.txt", params.indir.toString())}
                .combine(all_passed_id_channel)
                .filter{ meta, ani_best_hit, ids -> ids.contains(meta.id) }
                .map{ meta, ani_best_hit, ids -> [meta, ani_best_hit] }
                .groupTuple(by: 0)
                .map{ meta, files ->
                    def fileList = files.flatten()
                    if (fileList.size() > 1) {
                        log.info "--- [ANI DEBUG] --- Isolate: ${meta.id} | Found ${fileList.size()} files: ${fileList.collect { it.getName() }}"
                    }
                    def sortedFiles = fileList.sort { a, b ->
                        def dateA = (a.getName() =~ /\d{8}/) ? (a.getName() =~ /\d{8}/)[0].toInteger() : 0
                        def dateB = (b.getName() =~ /\d{8}/) ? (b.getName() =~ /\d{8}/)[0].toInteger() : 0
                        return dateB <=> dateA
                    }
                    def newest_file = sortedFiles[0]
                    if (fileList.size() > 1) {
                        log.info "--- [ANI DEBUG] --- Isolate: ${meta.id} | Selected Newest: ${newest_file.getName()}"
                    }
                    return [meta, newest_file]
                }

            // get prokka files
            def prokka_gff_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/annotation/*.gff')
            //create .gff and .faa files channel with meta information
            filtered_gff_ch = Channel.fromPath(prokka_gff_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta_non_extension(it,params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, gff, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, gff, all_passed_id_channel -> [meta, gff]} //remove all_passed_id_channel from output

            def prokka_faa_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/annotation/*.faa')
            filtered_faa_ch = Channel.fromPath(prokka_faa_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                .combine(all_passed_id_channel).filter{ meta, faa, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, faa, all_passed_id_channel -> [meta, faa]} //remove all_passed_id_channel from output

            def line_summary_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/*_summaryline.tsv')
            line_summary_ch = Channel.fromPath(line_summary_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_summaryline.tsv", params.indir.toString(),false)} // create meta for sample

            //////////////////////////// FOR CENTAR ///////////////////////////////////////

            // get files for MLST updating
            def combined_mlst_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/mlst/*_combined.tsv')

            //create .tax file channel with meta information
            filtered_combined_mlst_ch = Channel.fromPath(combined_mlst_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_combined.tsv", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, combined_mlst, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples - keep those in all_passed_id_channel
                .map{ meta, combined_mlst, all_passed_id_channel -> [meta, combined_mlst]} //remove all_passed_id_channel from output

            /////////////////////////// COLLECT SPECIES SPECIFIC FILES /////////////////////////////

            // get files for MLST updating
            def shigapass_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/ANI/*_ShigaPass_summary.csv')

            //collect .tax file channel with meta information
            shigapass_files_ch = Channel.fromPath(shigapass_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_ShigaPass_summary.csv", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, shiapass_files, all_passed_id -> all_passed_id.contains(meta.id)} //filtering out failured samples - keep those in all_passed_id_channel
                .map{ meta, shiapass_files, all_passed_id -> [meta, shiapass_files]} //remove all_passed_id_channel from output'


            //Will need to add in HV to the check eventually, but not today
            def gamma_hv_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/gamma_hv/*.gamma')
            //create .tax file channel with meta information
            filtered_gamma_hv_ch = Channel.fromPath(gamma_hv_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, gamma_hv, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples
                .map{ meta, gamma_hv, all_passed_id_channel -> [meta, gamma_hv]} //remove all_passed_id_channel from output

            // get .gamma files for AR updating
            def gamma_ar_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/gamma_ar/*ResGANNCBI*.gamma')

            if (params.mode_upper == "UPDATE_PHOENIX") {
                // 1. THE CANARY
                gamma_ar_with_flag = Channel.fromPath(gamma_ar_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, gamma_ar, ids -> ids.contains(meta.id)}
                    .map{ meta, gamma_ar, ids -> [meta, gamma_ar]}
                    .groupTuple()
                    .map{ meta, files -> [meta, files instanceof List ? files.flatten() : [files]] }
                    .combine(Channel.fromPath(params.ardb))
                    .map{ meta, gamma_ar, ardb -> check_update(meta, gamma_ar, ardb, "gamma", params.mode_upper) }
                    .multiMap { meta, file, flag ->
                        files: [meta, file]
                        flags: [meta, flag]
                    }

                filtered_gamma_ar_ch = gamma_ar_with_flag.files
                sample_needs_update_ch = gamma_ar_with_flag.flags

                // AMRFinder
                def armfinder_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/AMRFinder/*_all_genes{,_*}.tsv')
                filtered_amrfinder_ch = Channel.fromPath(armfinder_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, f, ids -> ids.contains(meta.id)}
                    .map{ meta, f, ids -> [meta, f]}
                    .groupTuple()
                    .map{ meta, files -> [meta, files instanceof List ? files.flatten() : [files]] }
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.amrfinder_db))
                    .map{ meta, f, needsUpdate, db -> check_update(meta, f, db, "amrfinder", params.mode_upper, needsUpdate) }
                    .map{ meta, file, flag -> [meta, file] }

                // Assembly Ratio
                def assembly_ratio_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/*_Assembly_ratio_*.txt')
                filtered_assembly_ratio_ch = Channel.fromPath(assembly_ratio_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, f, ids -> ids.contains(meta.id)}
                    .map{ meta, f, ids -> [meta, f]}
                    .groupTuple()
                    .map{ meta, files -> [meta, files instanceof List ? files.flatten() : [files]] }
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.ncbi_assembly_stats))
                    .map{ meta, f, needsUpdate, db -> check_update(meta, f, db, "ncbi_stats_ratio", params.mode_upper, needsUpdate) }
                    .map{ meta, file, flag -> [meta, file] }

                // GC Content
                def gc_content_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/*_GC_content_*.txt')
                filtered_gc_content_ch = Channel.fromPath(gc_content_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, f, ids -> ids.contains(meta.id)}
                    .map{ meta, f, ids -> [meta, f]}
                    .groupTuple()
                    .map{ meta, files -> [meta, files instanceof List ? files.flatten() : [files]] }
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.ncbi_assembly_stats))
                    .map{ meta, f, needsUpdate, db -> check_update(meta, f, db, "ncbi_stats_gc", params.mode_upper, needsUpdate) }
                    .map{ meta, file, flag -> [meta, file] }

                // SRST2
                def srst2_ar_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/srst2/*_srst2__results.txt')
                filtered_srst2_ar_ch = Channel.fromPath(srst2_ar_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, f, ids -> ids.contains(meta.id)}
                    .map{ meta, f, ids -> [meta, f]}
                    .groupTuple()
                    .map{ meta, files -> [meta, files instanceof List ? files.flatten() : [files]] }
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.ardb))
                    .map{ meta, f, needsUpdate, db -> check_update(meta, f, db, "srst2", params.mode_upper, needsUpdate) }
                    .map{ meta, file, flag -> [meta, file] }

                // Gamma PF
                def gamma_pf_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/gamma_pf/*PF-Replicons*.gamma')
                filtered_gamma_pf_ch = Channel.fromPath(gamma_pf_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, f, ids -> ids.contains(meta.id)}
                    .map{ meta, f, ids -> [meta, f]}
                    .groupTuple()
                    .map{ meta, files -> [meta, files instanceof List ? files.flatten() : [files]] }
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.gamdbpf))
                    .map{ meta, f, needsUpdate, db -> check_update(meta, f, db, "gamma_pf", params.mode_upper, needsUpdate) }
                    .map{ meta, file, flag -> [meta, file] }

            } else {
                // existing else branch - keep all the old logic unchanged
                filtered_gamma_ar_ch = Channel.fromPath(gamma_ar_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, gamma_ar, ids -> ids.contains(meta.id)}
                    .map{ meta, gamma_ar, ids -> [meta, gamma_ar]}
                    .combine(Channel.fromPath(params.ardb))
                    .map{ meta, gamma_ar, ardb ->
                        def ardbDate = ardb.getName() =~ /ResGANNCBI_(\d{8})_srst2\.fasta/
                        def matchingFile = gamma_ar.getName().contains(ardbDate[0][1]) ? gamma_ar : null
                        return [meta, matchingFile] }
                    .filter{ meta, gamma_ar -> gamma_ar != null }
                    .groupTuple()
                    .map { meta, files ->
                        def selected_file = files instanceof List ? files.sort().last() : files
                        return [ meta, selected_file ]
                    }

                def gamma_pf_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/gamma_pf/*PF-Replicons*.gamma')
                filtered_gamma_pf_ch = Channel.fromPath(gamma_pf_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, gamma_pf, ids -> ids.contains(meta.id)}
                    .map{ meta, gamma_pf, ids -> [meta, gamma_pf]}
                    .combine(Channel.fromPath(params.gamdbpf))
                    .map{ meta, gamma_pf, pfdb ->
                        def pfdbDate = pfdb.getName() =~ /PF-Replicons_(\d{8})\.fasta/
                        def matchingFile = gamma_pf.getName().contains(pfdbDate[0][1]) ? gamma_pf : null
                        return [meta, matchingFile] }
                    .filter{ meta, gamma_pf -> gamma_pf != null }

                def srst2_ar_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/srst2/*_srst2__results.txt')
                filtered_srst2_ar_ch = Channel.fromPath(srst2_ar_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, srst2_ar, ids -> ids.contains(meta.id)}
                    .map{ meta, srst2_ar, ids -> [meta, srst2_ar]}
                    .combine(Channel.fromPath(params.ardb))
                    .map{ meta, srst2_ar, ardb ->
                        def ardbDate = ardb.getName() =~ /ResGANNCBI_(\d{8})_srst2\.fasta/
                        def matchingFile = srst2_ar.getName().contains(ardbDate[0][1]) ? srst2_ar : null
                        return [meta, matchingFile] }
                    .filter{ meta, srst2_ar -> srst2_ar != null }

                def armfinder_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/AMRFinder/*_all_genes{,_*}.tsv')
                filtered_amrfinder_ch = Channel.fromPath(armfinder_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, f, ids -> ids.contains(meta.id)}
                    .map{ meta, f, ids -> [meta, f]}
                    .combine(Channel.fromPath(params.amrfinder_db))
                    .map{ meta, ncbi_report, ardb ->
                        def ardbDate = ardb.getName() =~ /amrfinderdb_v\d+\.\d+_(\d{8})\.\d+\.tar\.gz/
                        def matchingFile = ncbi_report.getName().contains(ardbDate[0][1]) ? ncbi_report : null
                        return [meta, matchingFile] }
                    .filter{ meta, ncbi_report -> ncbi_report != null }

                def assembly_ratio_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/*_Assembly_ratio_*.txt')
                filtered_assembly_ratio_ch = Channel.fromPath(assembly_ratio_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, ratio, ids -> ids.contains(meta.id)}
                    .map{ meta, ratio, ids -> [meta, ratio]}
                    .combine(Channel.fromPath(params.ncbi_assembly_stats))
                    .map{ meta, ratio, refdb ->
                        def refdbdate = (refdb.getName() =~ /_Assembly_stats_(\d{8})\.txt/)[0][1]
                        def matchingFile = ratio.getName().contains(refdbdate) ? ratio : null
                        return [meta, matchingFile] }
                    .filter{ meta, ratio -> ratio != null }
                    .groupTuple(by: 0)
                    .map{ meta, files ->
                        def newest_file = files.flatten().sort { a, b ->
                            def dateA = (a.getName() =~ /\d{8}/) ? (a.getName() =~ /\d{8}/)[0].toInteger() : 0
                            def dateB = (b.getName() =~ /\d{8}/) ? (b.getName() =~ /\d{8}/)[0].toInteger() : 0
                            return dateB <=> dateA
                        }[0]
                        return [meta, newest_file]
                    }

                def gc_content_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/*_GC_content_*.txt')
                filtered_gc_content_ch = Channel.fromPath(gc_content_glob)
                    .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())}
                    .combine(all_passed_id_channel).filter{ meta, gc_content, ids -> ids.contains(meta.id)}
                    .map{ meta, gc_content, ids -> [meta, gc_content]}
                    .combine(Channel.fromPath(params.ncbi_assembly_stats))
                    .map{ meta, gc_content, refdb ->
                        def refdbdate = refdb.getName() =~ /_Assembly_stats_(\d{8})\.txt/
                        def matchingFile = gc_content.getName().contains(refdbdate[0][1]) ? gc_content : null
                        return [meta, matchingFile] }
                    .filter{ meta, gc_content -> gc_content != null }
            }


            // get CENTAR files
            def centar_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/CENTAR/*_centar_output.tsv')

            //collect .tax file channel with meta information
            centar_files_ch = Channel.fromPath(centar_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_centar_output.tsv", params.indir.toString(),false)} // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, centar_files, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)} //filtering out failured samples - keep those in all_passed_id_channel
                .map{ meta, centar_files, all_passed_id_channel -> [meta, centar_files]}.ifEmpty( [[id: "", project_id: ""], []] ) //remove all_passed_id_channel from output

            /////////////////////////// COLLECT README FOR UPDATER ///////////////////////////////
            def readme_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/*_updater_log.tsv')

            readme_files_ch = Channel.fromPath(readme_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta(it, "_updater_log.tsv", params.indir.toString(), false)}.ifEmpty( [[id: "", project_id: ""], []] )  // create meta for sample
                .combine(all_passed_id_channel).filter{ meta, readme_files, all_passed_id_channel -> all_passed_id_channel.contains(meta.id)}.ifEmpty( [[id: "", project_id: ""], [], []] ) //filtering out failed samples - keep those in all_passed_id_channel
                .map{ meta, readme_files, all_passed_id_channel -> [meta, readme_files]} //remove all_passed_id_channel from output

            /////////////////////////// COLLECT PROJECT LEVEL FILES ///////////////////////////////

            def synopsis_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*/*.synopsis')
            synopsis_ch = Channel.fromPath(synopsis_glob) // use created regrex to get samples
                .map{ it -> InputChannelUtils.create_meta_non_extension(it, params.indir.toString())} // create meta for sample and adding group ID
            //filtering out failured samples
            filtered_synopsis_ch = synopsis_ch.filter{ meta, synopsis -> all_passed_id_channel.contains(meta.id)}

            // use created regrex to get samples
            def griphin_excel_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*_GRiPHin_Summary.xlsx')
            all_griphin_excel_ch = all_passed_id_channel.flatten().combine(Channel.fromPath(griphin_excel_glob)).map{ it -> InputChannelUtils.create_groups_and_id(it, params.indir.toString())}
            //create regrex, get files in dir, add in meta information, change name
            def griphin_tsv_glob = InputChannelUtils.append_to_path(params.indir.toString(),'*_GRiPHin_Summary.tsv')
            all_griphin_tsv_ch = all_passed_id_channel.flatten().combine(Channel.fromPath(griphin_tsv_glob)).map{ it -> InputChannelUtils.create_groups_and_id(it, params.indir.toString())}
            def phoenix_tsv_glob = InputChannelUtils.append_to_path(params.indir.toString(),'Phoenix_Summary.tsv')
            all_phoenix_tsv_ch = all_passed_id_channel.flatten().combine(Channel.fromPath(phoenix_tsv_glob)).map{ it -> InputChannelUtils.create_groups_id_and_busco(it, params.indir.toString())}
            def pipeline_info_glob = InputChannelUtils.append_to_path(params.indir.toString(),'pipeline_info/software_versions.yml')
            // Keep the original logic for the Project-level tools
            all_pipeline_info_ch = all_passed_id_channel.flatten()
                .combine(Channel.fromPath(pipeline_info_glob))
                .map{ it -> InputChannelUtils.create_groups_and_id(it, params.indir.toString())}

            // 1. First, define where the update info should be
            def update_info_path = "${params.indir}/update_pipeline_info/software_versions.yml"
            def update_file = file(update_info_path)

            // 2. Explicitly create the channel so it EXISTS for the join
            // If the file is there, use it; if not, send an empty list
            update_pipeline_info_ch = all_passed_id_channel
                .map { id -> 
                    def meta = [project_id: params.indir.toString().split('/')[-1]]
                    return [ meta, update_file.exists() ? update_file : [] ]
                }
                .unique()

            CREATE_SAMPLESHEET (
                indir
            )
            ch_versions = ch_versions.mix(CREATE_SAMPLESHEET.out.versions)

            valid_samplesheet = CREATE_SAMPLESHEET.out.samplesheet

            summary_files_ch = valid_samplesheet.flatten().splitCsv( header:true, sep:',' )
                .map{ it -> InputChannelUtils.create_summary_files_channels(it) }
                .map{ meta, griphin_excel, griphin_tsv, phoenix_tsv, pipeline_info, update_info -> 
                    [[project_id: meta.project_id], griphin_excel, griphin_tsv, phoenix_tsv, pipeline_info, update_info ?: []]
                }
                .unique()
                .collect()

            // combining all summary files into one channel
            summary_files_ch = all_griphin_excel_ch.map{            meta, griphin_excel -> [[project_id:meta.project_id], griphin_excel]}
                        .join(all_griphin_tsv_ch.map{   meta, griphin_tsv   -> [[project_id:meta.project_id], griphin_tsv]},   by: [0])
                        .join(all_phoenix_tsv_ch.map{   meta, phoenix_tsv   -> [[project_id:meta.project_id], phoenix_tsv]},   by: [0])
                        .join(all_pipeline_info_ch.map{ meta, pipeline_info -> [[project_id:meta.project_id], pipeline_info]}, by: [0])
                        .join(update_pipeline_info_ch.map{ meta, update_info -> [[project_id:meta.project_id], update_info]}, by: [0], remainder: true)
                        .map{ meta, griphin_excel, griphin_tsv, phoenix_tsv, pipeline_info, update_info -> 
                            [[project_id:meta.project_id.toString().split('/')[-1].replace("]", "")], griphin_excel, griphin_tsv, phoenix_tsv, pipeline_info, update_info ?: []] 
                        }
                        .view { "BEFORE UNIQUE: ${it[0].project_id} | Excel:${it[1]?.name} | TSV:${it[2]?.name} | Phoenix:${it[3]?.name} | Info:${it[4]?.name} | Update:${it[5]?.name}" }
                        .unique { it[0].project_id }
                        .view { "AFTER UNIQUE: ${it[0].project_id} | Excel:${it[1]?.name} | TSV:${it[2]?.name} | Phoenix:${it[3]?.name} | Info:${it[4]?.name} | Update:${it[5]?.name}" }

            summary_files_ch.view { meta, excel, g_tsv, p_tsv, p_info, u_info ->
                """
                DEBUG summary_files_ch:
                Project ID: ${meta.project_id}
                Excel:      ${excel.name}
                Griphin TSV:${g_tsv.name}
                Phoenix TSV:${p_tsv.name}
                P-Info:     ${p_info.name}
                U-Info:     ${u_info ? u_info.name : 'EMPTY'}
                -------------------------------------------
                """
            }

            // pulling all the necessary project level files into channels
            COLLECT_PROJECT_FILES (
                summary_files_ch, Channel.value(false)
            )
            ch_versions = ch_versions.mix(COLLECT_PROJECT_FILES.out.versions)

            // 2. Assign the project-level variables FIRST
            pipeline_info_ch        = COLLECT_PROJECT_FILES.out.software_versions_file
            update_pipeline_info_ch = COLLECT_PROJECT_FILES.out.update_software_versions

            // Correcting the structure so it's [ [meta], file ] 
            // instead of [ id, [meta], file ]
            // Fix for Closure 163 (Pipeline Info)
            isolate_version_broadcast_ch = all_passed_id_channel.flatten()
                .combine(pipeline_info_ch) 
                .map { id, project_meta, file -> 
                    // Merge the ID into the project meta map
                    def new_meta = project_meta + [id: id]
                    return [ new_meta, file ] 
                }

            isolate_update_broadcast_ch = all_passed_id_channel.flatten()
                .combine(update_pipeline_info_ch)
                .map { id, project_meta, file -> 
                    def new_meta = project_meta + [id: id]
                    return [ new_meta, file ] 
                }

            griphin_excel_ch = COLLECT_PROJECT_FILES.out.griphin_excel
            griphin_tsv_ch = COLLECT_PROJECT_FILES.out.griphin_tsv
            phoenix_tsv_ch = COLLECT_PROJECT_FILES.out.phoenix_tsv.map{it -> InputChannelUtils.add_entry_meta(it)}
            pipeline_info_ch = COLLECT_PROJECT_FILES.out.software_versions_file
            sample_pipeline_info_ch = COLLECT_PROJECT_FILES.out.software_versions_file
            update_pipeline_info_ch = COLLECT_PROJECT_FILES.out.update_software_versions
            samplesheet_meta_ch = Channel.empty().ifEmpty([]) // No multiple samplesheets meta channel for indir

        } else if (samplesheet != null) {

            def ensureSingle = { meta, files ->
                if (files instanceof List) {
                    if (files.size() == 0) return [meta, []]
                    def single = files.flatten().sort { a, b ->
                        def dA = (a.name =~ /\d{8}/) ? (a.name =~ /\d{8}/)[0].toInteger() : 0
                        def dB = (b.name =~ /\d{8}/) ? (b.name =~ /\d{8}/)[0].toInteger() : 0
                        return dB <=> dA
                    }[0]
                    return [meta, single]
                }
                return [meta, files]
            }

            meta_ch = Channel.fromPath(samplesheet).splitCsv( header:true, sep:',' ).map{ InputChannelUtils.create_samplesheet_meta(it) }.unique()

            SAMPLESHEET_CHECK (
                samplesheet, false, false, true, false, meta_ch
            )
            ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

            samplesheet = SAMPLESHEET_CHECK.out.csv.first()
            samplesheet_meta_ch = SAMPLESHEET_CHECK.out.csv_by_dir.flatten().map{ it -> InputChannelUtils.transformSamplesheets(it)}

            // To make things backwards compatible we need to check if the file_integrity sample is there and if not create it.
            file_integrity_exists = samplesheet.splitCsv( header:true, sep:',' ).map{ it -> InputChannelUtils.check_file_integrity(it) }.filter{meta, clean_path, fairy_exists -> fairy_exists == false }

            CREATE_FAIRY_FILE (
                file_integrity_exists, false
            )
            ch_versions = ch_versions.mix(CREATE_FAIRY_FILE.out.versions)

            directory_ch = samplesheet.splitCsv( header:true, sep:',' ).map{ it -> InputChannelUtils.create_dir_channels(it) }

            //adding meta.id to end of dir - otherwise too many files are copied and it takes forever.
            sample_directory_ch = samplesheet.splitCsv( header:true, sep:',' ).map{ it -> InputChannelUtils.create_sample_dir_channels(it) }

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
            filtered_renamed_scaffolds_ch = COLLECT_SAMPLE_FILES.out.renamed_scaffolds
            filtered_scaffolds_ch         = COLLECT_SAMPLE_FILES.out.scaffolds
            filtered_gff_ch               = COLLECT_SAMPLE_FILES.out.gff
            filtered_faa_ch               = COLLECT_SAMPLE_FILES.out.faa
            line_summary_ch               = COLLECT_SAMPLE_FILES.out.summary_line
            filtered_synopsis_ch          = COLLECT_SAMPLE_FILES.out.synopsis
            filtered_taxonomy_ch          = COLLECT_SAMPLE_FILES.out.tax
            filtered_gamma_hv_ch          = COLLECT_SAMPLE_FILES.out.gamma_hv
            if (params.mode_upper == "UPDATE_PHOENIX") {
                // 1. THE CANARY: Run the check for Gamma
                gamma_ar_with_flag = COLLECT_SAMPLE_FILES.out.gamma_ar
                    .combine(Channel.fromPath(params.ardb))
                    .map{ meta, gamma_ar, ardb -> 
                        check_update(meta, gamma_ar, ardb, "gamma", params.mode_upper) 
                    }
                    .multiMap { meta, file, flag ->
                        files: [meta, file]
                        flags: [meta, flag]
                    }

                filtered_gamma_ar_ch = gamma_ar_with_flag.files

                sample_needs_update_ch = gamma_ar_with_flag.flags
                    .map{ meta, flag ->
//                        System.err.println ">>> FLAG KEY [${meta.id}]: project_id=${meta.project_id} flag=${flag}"
                        [meta, flag]
                    }

/*                COLLECT_SAMPLE_FILES.out.amrfinder_report.map { meta, f ->
                    log.info ">>> AMRFINDER_REPORT: id=${meta.id} | file=${f}"
                    [meta, f]
                }

                sample_needs_update_ch.map { meta, flag ->
                    log.info ">>> SAMPLE_NEEDS_UPDATE: id=${meta.id} | flag=${flag}"
                    [meta, flag]
                }
*/
                 // AMRFinder
                filtered_amrfinder_ch = COLLECT_SAMPLE_FILES.out.amrfinder_report
                    .map{ meta, f -> 
//                        System.err.println ">>> AMRFINDER KEY [${meta.id}]: project_id=${meta.project_id}"
                        [meta, f] 
                    }
                    .join(sample_needs_update_ch) 
                    .combine(Channel.fromPath(params.amrfinder_db))
                    .map{ meta, f, needsUpdate, db -> 
                        check_update(meta, f, db, "amrfinder", params.mode_upper, needsUpdate) 
                    }
                    .map{ meta, file, flag -> [meta, file] }

                // Assembly Ratio
                filtered_assembly_ratio_ch = COLLECT_SAMPLE_FILES.out.assembly_ratio
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.ncbi_assembly_stats))
                    .map{ meta, f, needsUpdate, db -> 
                        check_update(meta, f, db, "ncbi_stats_ratio", params.mode_upper, needsUpdate) 
                    }
                    .map{ meta, file, flag -> [meta, file] }

                // GC Content
                filtered_gc_content_ch = COLLECT_SAMPLE_FILES.out.gc_content
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.ncbi_assembly_stats))
                    .map{ meta, f, needsUpdate, db -> 
                        check_update(meta, f, db, "ncbi_stats_gc", params.mode_upper, needsUpdate) 
                    }
                    .map{ meta, file, flag -> [meta, file] }

                // SRST2
                filtered_srst2_ar_ch = COLLECT_SAMPLE_FILES.out.srst2_ar
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.ardb))
                    .map{ meta, f, needsUpdate, db -> 
                        check_update(meta, f, db, "srst2", params.mode_upper, needsUpdate) 
                    }
                    .map{ meta, file, flag -> [meta, file] }

                // Gamma PF
                filtered_gamma_pf_ch = COLLECT_SAMPLE_FILES.out.gamma_pf
                    .join(sample_needs_update_ch)
                    .combine(Channel.fromPath(params.gamdbpf))
                    .map{ meta, f, needsUpdate, db -> 
                        check_update(meta, f, db, "gamma_pf", params.mode_upper, needsUpdate) 
                    }
                    .map{ meta, file, flag -> [meta, file] }          

            } else {
                filtered_gamma_ar_ch       = COLLECT_SAMPLE_FILES.out.gamma_ar.map{ m, f -> ensureSingle(m, f) }
                filtered_amrfinder_ch      = COLLECT_SAMPLE_FILES.out.amrfinder_report.map{ m, f -> ensureSingle(m, f) }
                filtered_assembly_ratio_ch = COLLECT_SAMPLE_FILES.out.assembly_ratio.map{ m, f -> ensureSingle(m, f) }
                filtered_gc_content_ch     = COLLECT_SAMPLE_FILES.out.gc_content.map{ m, f -> ensureSingle(m, f) }
                filtered_srst2_ar_ch       = COLLECT_SAMPLE_FILES.out.srst2_ar.map{ m, f -> ensureSingle(m, f) }
                filtered_gamma_pf_ch       = COLLECT_SAMPLE_FILES.out.gamma_pf.map{ m, f -> ensureSingle(m, f) }
            }
            filtered_trimd_kraken_bh_ch        = COLLECT_SAMPLE_FILES.out.trimd_kraken_bh
            filtered_trimd_krona_ch            = COLLECT_SAMPLE_FILES.out.trimd_kraken_krona
            filtered_trimd_kraken_report_ch    = COLLECT_SAMPLE_FILES.out.trimd_kraken_report
            filtered_wtasmbld_kraken_bh_ch     = COLLECT_SAMPLE_FILES.out.wtasmbld_kraken_bh
            filtered_wtasmbld_krona_ch         = COLLECT_SAMPLE_FILES.out.wtasmbld_kraken_krona
            filtered_wtasmbld_kraken_report_ch = COLLECT_SAMPLE_FILES.out.wtasmbld_kraken_report
            filtered_asmbld_kraken_bh_ch       = COLLECT_SAMPLE_FILES.out.asmbld_kraken_bh
            filtered_asmbld_krona_ch           = COLLECT_SAMPLE_FILES.out.asmbld_kraken_krona
            filtered_asmbld_kraken_report_ch   = COLLECT_SAMPLE_FILES.out.asmbld_kraken_report
            filtered_busco_short_summary_ch    = COLLECT_SAMPLE_FILES.out.busco_short_summary
            filtered_trimmed_stats_ch          = COLLECT_SAMPLE_FILES.out.trimmed_stats
            filtered_raw_stats_ch              = COLLECT_SAMPLE_FILES.out.raw_stats
            filtered_quast_ch                  = COLLECT_SAMPLE_FILES.out.quast_report
            filtered_ani_ch                    = COLLECT_SAMPLE_FILES.out.ani
            filtered_ani_best_hit_ch           = COLLECT_SAMPLE_FILES.out.ani_best_hit.map{ m, f -> ensureSingle(m, f) }
            filtered_combined_mlst_ch          = COLLECT_SAMPLE_FILES.out.combined_mlst

            //species specific files
            shigapass_files_ch = COLLECT_SAMPLE_FILES.out.shigapass_output
            centar_files_ch = COLLECT_SAMPLE_FILES.out.centar_output
            //readme files
            readme_files_ch = COLLECT_SAMPLE_FILES.out.readme

            summary_files_ch = samplesheet.flatten().splitCsv( header:true, sep:',' )
                .map{ it -> InputChannelUtils.create_summary_files_channels(it) }
                .map{ meta, griphin_excel, griphin_tsv, phoenix_tsv, pipeline_info, update_info -> 
                    [[project_id: meta.project_id], griphin_excel, griphin_tsv, phoenix_tsv, pipeline_info, update_info ?: []]
                }
                .unique()

            COLLECT_PROJECT_FILES (
                summary_files_ch, Channel.value(true)
            )
            ch_versions = ch_versions.mix(COLLECT_PROJECT_FILES.out.versions)

            griphin_excel_ch = COLLECT_PROJECT_FILES.out.griphin_excel
            griphin_tsv_ch = COLLECT_PROJECT_FILES.out.griphin_tsv
            phoenix_tsv_ch = COLLECT_PROJECT_FILES.out.phoenix_tsv.map{it -> InputChannelUtils.add_entry_meta(it)}
            pipeline_info_ch = COLLECT_PROJECT_FILES.out.software_versions_file
            update_pipeline_info_ch = COLLECT_PROJECT_FILES.out.update_software_versions

            valid_samplesheet = samplesheet

            isolate_version_broadcast_ch = filtered_scaffolds_ch
                .map { meta, files -> [[project_id: meta.project_id], meta] }
                .combine(pipeline_info_ch.map { meta, file -> [[project_id: meta.project_id], file] }, by: 0)
                .map { project_key, sample_meta, pipeline_info -> [sample_meta, pipeline_info] }

            sample_pipeline_info_ch = filtered_scaffolds_ch
                .map { meta, files -> [[project_id: meta.project_id], meta] }
                .combine(pipeline_info_ch.map { meta, file -> [[project_id: meta.project_id], file] }, by: 0)
                .map { project_meta, sample_meta, pipeline_info -> [sample_meta, pipeline_info] }

            isolate_update_broadcast_ch = filtered_scaffolds_ch
                .map { meta, files -> [[project_id: meta.project_id], meta] }
                .combine(update_pipeline_info_ch.map { meta, file -> [[project_id: meta.project_id], file] }, by: 0)
                .map { project_key, sample_meta, update_pipeline_info -> [sample_meta, update_pipeline_info] }

            sample_update_pipeline_info_ch = filtered_scaffolds_ch
                .map { meta, files -> [[project_id: meta.project_id], meta] }
                .combine(update_pipeline_info_ch.map { meta, file -> [[project_id: meta.project_id], file] }, by: 0)
                .map { project_meta, sample_meta, update_pipeline_info -> [sample_meta, update_pipeline_info] }

        } else {
            exit 1, 'You need EITHER an input directory samplesheet (using --input) or a single sample directory (using --indir)!'
        }

    emit:
        //project level summary files
        griphin_excel_ch      = griphin_excel_ch
        griphin_tsv_ch        = griphin_tsv_ch
        phoenix_tsv_ch        = phoenix_tsv_ch
        pipeline_info_ch      = pipeline_info_ch
        directory_ch          = directory_ch
        valid_samplesheet     = valid_samplesheet
        versions              = ch_versions
        pipeline_info         = pipeline_info_ch
        update_pipeline_info  = update_pipeline_info_ch


        //species specific files
        centar             = centar_files_ch
        shigapass          = shigapass_files_ch
        //updater
        readme             = readme_files_ch
        samplesheet_meta_ch = samplesheet_meta_ch

        // sample specific files
        renamed_scaffolds      = filtered_renamed_scaffolds_ch
        filtered_scaffolds     = filtered_scaffolds_ch      // channel: [ meta, [ scaffolds_file ] ]
        reads                  = combined_reads_ch
        taxonomy               = filtered_taxonomy_ch
        prokka_gff             = filtered_gff_ch
        prokka_faa             = filtered_faa_ch
        fairy_outcome          = file_integrity_ch
        line_summary           = line_summary_ch // need non-filtered to make summary files will all samples in project folder
        synopsis               = filtered_synopsis_ch
        ani                    = filtered_ani_ch
        ani_best_hit           = filtered_ani_best_hit_ch
        ncbi_report            = filtered_amrfinder_ch
        gamma_ar               = filtered_gamma_ar_ch
        srst2_ar               = filtered_srst2_ar_ch
        gamma_pf               = filtered_gamma_pf_ch
        gamma_hv               = filtered_gamma_hv_ch
        assembly_ratio         = filtered_assembly_ratio_ch
        gc_content             = filtered_gc_content_ch
        k2_trimd_bh_summary    = filtered_trimd_kraken_bh_ch
        k2_trimd_krona         = filtered_trimd_krona_ch
        k2_trimd_report        = filtered_trimd_kraken_report_ch
        k2_wtasmbld_bh_summary = filtered_wtasmbld_kraken_bh_ch
        k2_wtasmbld_krona      = filtered_wtasmbld_krona_ch
        k2_wtasmbld_report     = filtered_wtasmbld_kraken_report_ch
        k2_asmbld_bh_summary   = filtered_asmbld_kraken_bh_ch
        k2_asmbld_krona        = filtered_asmbld_krona_ch
        k2_asmbld_report       = filtered_asmbld_kraken_report_ch
        busco_short_summary    = filtered_busco_short_summary_ch
        fastp_total_qc         = filtered_trimmed_stats_ch
        raw_stats              = filtered_raw_stats_ch
        quast_report           = filtered_quast_ch
        combined_mlst          = filtered_combined_mlst_ch // for centar entry
        pipeline_info_isolate  = isolate_version_broadcast_ch // Use this for summary lines
        update_pipeline_info_isolate = isolate_update_broadcast_ch
        sample_needs_update_ch  // already exists, just needs to be added to emit block
}
