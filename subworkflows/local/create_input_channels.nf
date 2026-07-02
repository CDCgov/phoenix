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

def newest_by_embedded_date = { ch ->
    ch.groupTuple(by: 0)
        .map { meta, files ->
            def newest = files.flatten().sort { a, b ->
                def da = (a.getName() =~ /\d{8}/) ? (a.getName() =~ /\d{8}/)[0].toInteger() : 0
                def db = (b.getName() =~ /\d{8}/) ? (b.getName() =~ /\d{8}/)[0].toInteger() : 0
                db <=> da
            }[0]
            [meta, newest]
        }
}



workflow CREATE_INPUT_CHANNELS {
    take:
        samplesheet  // params.input
        centar       // true when centar is run

    main:

        ch_versions = Channel.empty() // Used to collect the software versions

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
        //file_integrity_ch = CREATE_FAIRY_FILE.out.created_fairy_file.collect().ifEmpty([]).combine(COLLECT_SAMPLE_FILES.out.fairy_summary.collect().ifEmpty([])).flatten().buffer(size:2)
        file_integrity_ch = COLLECT_SAMPLE_FILES.out.fairy_summary
            .concat(
                CREATE_FAIRY_FILE.out.created_fairy_file
            )
            .unique { meta, files -> meta.id }

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
            fairy_passed_scaffolds_ch = COLLECT_SAMPLE_FILES.out.scaffolds
                .join(file_integrity_ch, by: [0])
                .filter { meta, scaffolds, fairy_outcome ->
                    def target_file = fairy_outcome instanceof List 
                        ? fairy_outcome.find { it.name.endsWith('_scaffolds_summary.txt') } 
                        : fairy_outcome
                    def content = file(target_file.toString()).text
                    !content.contains("FAILED")
                }
                .map { meta, scaffolds, fairy_outcome -> [meta, scaffolds] }
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

            // Find samples with scaffolds but NO gamma_ar file - these always need updating
            missing_gamma_ar_ch = filtered_scaffolds_ch
                .map{ meta, scaffolds -> [meta.id, meta] }
                .join(
                    gamma_ar_with_flag.flags.map{ meta, flag -> [meta.id, flag] },
                    remainder: true
                )
                .filter{ id, meta, flag -> flag == null }
                .map{ id, meta, flag -> [meta, true] }

            sample_needs_update_ch = gamma_ar_with_flag.flags
                .mix(missing_gamma_ar_ch)

            // GAMMA HV
            filtered_gamma_hv_ch = COLLECT_SAMPLE_FILES.out.gamma_hv
                .map{ meta, f -> 
                    [meta, f] 
                }
                .join(sample_needs_update_ch) 
                .combine(Channel.fromPath(params.hvgamdb))
                .map{ meta, f, needsUpdate, db -> 
                    check_update(meta, f, db, "gamma_hv", params.mode_upper, needsUpdate) 
                }
                .map{ meta, file, flag -> [meta, file] }

            // AMRFinder
            filtered_amrfinder_ch = COLLECT_SAMPLE_FILES.out.amrfinder_report
                .map{ meta, f -> 
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
            sample_needs_update_ch     = Channel.empty() // emit empty channel if not update mode
            fairy_passed_scaffolds_ch  = Channel.empty() // emit empty channel if not update mode
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
        filtered_ani_ch                    = COLLECT_SAMPLE_FILES.out.ani.map{ m, f -> ensureSingle(m, f) }
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

        // Function to detect what the original run type was based on what files are present in the directory
        mode_type_ch = directory_ch.map { meta, dir ->
            [meta, InputChannelUtils.detect_mode_type(file("${dir}/${meta.id}"))]
        }

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
        sample_needs_update_ch
        mode_type = mode_type_ch
        fairy_passed_scaffolds = fairy_passed_scaffolds_ch
}
