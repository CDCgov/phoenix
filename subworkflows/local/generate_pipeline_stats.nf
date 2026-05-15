//
// Subworkflow: Running SPAdes and checking if spades failed to create scaffolds
//

include { GENERATE_PIPELINE_STATS } from '../../modules/local/generate_pipeline_stats'

// Groovy funtion to make [ meta.id, [] ] - just an empty channel
def create_empty_ch(input_for_meta) { // We need meta.id associated with the empty list which is why .ifempty([]) won't work
    def meta_id
    meta_id = input_for_meta[0]
    def output_array
    output_array = [ meta_id, [] ]
    return output_array
}

workflow GENERATE_PIPELINE_STATS_WF {
    take:
        fastp_raw_qc           // channel: tuple (meta) path(fastp_raw_qc): GATHERING_READ_QC_STATS.out.fastp_raw_qc
        fastp_total_qc         // channel: tuple (meta) path(fastp_total_qc): GATHERING_READ_QC_STATS.out.fastp_total_qc
        fullgene_results       // channel: tuple (meta) path(fullgene_results): SRST2_TRIMD_AR.out.fullgene_results
        trimd_report           // channel: tuple (meta) path(report): KRAKEN2_TRIMD.out.report
        trimd_krona_html       // channel: tuple (meta) path(krona_html): KRAKEN2_TRIMD.out.krona_html
        trimd_k2_bh_summary    // channel: tuple (meta) path(k2_bh_summary): KRAKEN2_TRIMD.out.k2_bh_summary
        renamed_fastas
        filtered_fastas
        mlst
        gamma_hv
        gamma_ar
        gamma_pf
        quast_report
        busco
        asmbld_report          // channel: tuple (meta) path(report): KRAKEN2_ASMBLD.out.report
        asmbld_krona_html      // channel: tuple (meta) path(krona_html): KRAKEN2_ASMBLD.out.krona_html
        asmbld_k2_bh_summary   // channel: tuple (meta) path(k2_bh_summary): KRAKEN2_ASMBLD.out.k2_bh_summary
        wtasmbld_report        // channel: tuple (meta) path(report): KRAKEN2_WTASMBLD.out.report
        wtasmbld_krona_html    // channel: tuple (meta) path(krona_html): KRAKEN2_WTASMBLD.out.krona_html
        wtasmbld_k2_bh_summary // channel: tuple (meta) path(k2_bh_summary): KRAKEN2_WTASMBLD.out.k2_bh_summary
        taxa_id
        format_ani
        assembly_ratio
        amr_point_mutations    // channel: tuple val(meta), path(report): AMRFINDERPLUS_RUN.out.report
        gc_content             // CALCULATE_ASSEMBLY_RATIO.out.gc_content
        run_type               // Easy enough for standard modes. Updater will provide the preivous run_type through the create_input_channels variable
        //in_extended_qc       // true for internal phoenix and false otherwise

    main:
        ch_versions = Channel.empty() // Used to collect the software versions
        def pipeline_stats_ch

        // Null-guard all inputs — replace null with empty channel
        if (wtasmbld_report == null)      wtasmbld_report      = Channel.empty()
        if (run_type == null)             run_type             = Channel.empty()
        if (fastp_raw_qc == null)         fastp_raw_qc         = Channel.empty()
        if (renamed_fastas == null)       renamed_fastas       = Channel.empty()
        if (filtered_fastas == null)      filtered_fastas      = Channel.empty()
        if (mlst == null)                 mlst                 = Channel.empty()
        if (busco == null)                busco                = Channel.empty()
        if (asmbld_report == null)        asmbld_report        = Channel.empty()
        if (asmbld_krona_html == null)    asmbld_krona_html    = Channel.empty()
        if (asmbld_k2_bh_summary == null) asmbld_k2_bh_summary = Channel.empty()
        if (fullgene_results == null)     fullgene_results     = Channel.empty()

        def add_padding = { ch, id_ch ->
            ch.mix(
                wtasmbld_report
                    .map { meta, report -> tuple(meta.id, meta, report) }
                    .join(id_ch, by: 0)
                    .map { id, meta, report, flag -> [meta, []] }
            )
            .unique { row -> [row[0].id, row[0].project_id] }
        }

        wtasmbld_report_with_rt = wtasmbld_report
            .join(run_type, by: [0])  // emits [meta, report, rt]

        // Collect sample IDs that need fullgene/SRST2 padding (everything except CDC_PHOENIX)
        no_fullgene_ids = wtasmbld_report_with_rt
            .filter { meta, report, rt ->
                rt.base == "PHOENIX" || rt.base == "SCAFFOLDS" || rt.base == "CDC_SCAFFOLDS"
            }
            .map { meta, report, rt -> tuple(meta.id, true) }

        // Collect sample IDs that need asmbld/busco padding (PHOENIX and SCAFFOLDS only)
        no_asmbld_kraken_ids = wtasmbld_report_with_rt
            .filter { meta, report, rt ->
                rt.base == "PHOENIX" || rt.base == "SCAFFOLDS"
            }
            .map { meta, report, rt -> tuple(meta.id, true) }

        // Collect scaffold-origin sample IDs that need reads/trimd padding (SCAFFOLDS and CDC_SCAFFOLDS)
        scaffold_origin_ids = wtasmbld_report_with_rt
            .filter { meta, report, rt ->
                rt.base == "SCAFFOLDS" || rt.base == "CDC_SCAFFOLDS"
            }
            .map { meta, report, rt -> tuple(meta.id, true) }

        // Strip run_type back off after identifying scaffold-origin samples
        wtasmbld_report = wtasmbld_report_with_rt
            .map { meta, report, rt ->
                [meta, report]
            }

        if (params.mode_upper == "SCAFFOLDS") {
            // All samples are scaffold-based — pad all trimmed read channels wholesale
            fastp_raw_qc        = wtasmbld_report.map{ it -> create_empty_ch(it) }
            fastp_total_qc      = wtasmbld_report.map{ it -> create_empty_ch(it) }
            fullgene_results    = wtasmbld_report.map{ it -> create_empty_ch(it) }
            trimd_report        = wtasmbld_report.map{ it -> create_empty_ch(it) }
            trimd_krona_html    = wtasmbld_report.map{ it -> create_empty_ch(it) }
            trimd_k2_bh_summary = wtasmbld_report.map{ it -> create_empty_ch(it) }
            // Add these:
            busco                = wtasmbld_report.map{ it -> create_empty_ch(it) }
            asmbld_report        = wtasmbld_report.map{ it -> create_empty_ch(it) }
            asmbld_krona_html    = wtasmbld_report.map{ it -> create_empty_ch(it) }
            asmbld_k2_bh_summary = wtasmbld_report.map{ it -> create_empty_ch(it) }
        }

        if (params.mode_upper == "CDC_SCAFFOLDS") {
            // All samples are scaffold-based — pad all trimmed read channels wholesale
            fastp_raw_qc        = wtasmbld_report.map{ it -> create_empty_ch(it) }
            fastp_total_qc      = wtasmbld_report.map{ it -> create_empty_ch(it) }
            fullgene_results    = wtasmbld_report.map{ it -> create_empty_ch(it) }
            trimd_report        = wtasmbld_report.map{ it -> create_empty_ch(it) }
            trimd_krona_html    = wtasmbld_report.map{ it -> create_empty_ch(it) }
            trimd_k2_bh_summary = wtasmbld_report.map{ it -> create_empty_ch(it) }
        }

        if (params.mode_upper == "PHOENIX" ) {
            busco            = wtasmbld_report.map{ it -> create_empty_ch(it) }
            asmbld_report    = wtasmbld_report.map{ it -> create_empty_ch(it) }
            asmbld_krona_html    = wtasmbld_report.map{ it -> create_empty_ch(it) }
            asmbld_k2_bh_summary = wtasmbld_report.map{ it -> create_empty_ch(it) }
            fullgene_results = wtasmbld_report.map{ it -> create_empty_ch(it) }
        }

        if (params.mode_upper == "UPDATE_PHOENIX") {
            fastp_raw_qc         = add_padding(fastp_raw_qc,         scaffold_origin_ids)
            fastp_total_qc       = add_padding(fastp_total_qc,        scaffold_origin_ids)
            trimd_report         = add_padding(trimd_report,          scaffold_origin_ids)
            trimd_krona_html     = add_padding(trimd_krona_html,      scaffold_origin_ids)
            trimd_k2_bh_summary  = add_padding(trimd_k2_bh_summary,  scaffold_origin_ids)
            fullgene_results     = add_padding(fullgene_results,      no_fullgene_ids)
            busco                = add_padding(busco,                 no_asmbld_kraken_ids)
            asmbld_report        = add_padding(asmbld_report,         no_asmbld_kraken_ids)
            asmbld_krona_html    = add_padding(asmbld_krona_html,     no_asmbld_kraken_ids)
            asmbld_k2_bh_summary = add_padding(asmbld_k2_bh_summary, no_asmbld_kraken_ids)
            mlst                 = add_padding(mlst,                  scaffold_origin_ids)
            gamma_hv             = add_padding(gamma_hv,              scaffold_origin_ids)
            gamma_ar             = add_padding(gamma_ar,              scaffold_origin_ids)
            gamma_pf             = add_padding(gamma_pf,              scaffold_origin_ids)
            renamed_fastas       = add_padding(renamed_fastas,        scaffold_origin_ids)
            filtered_fastas      = add_padding(filtered_fastas,       scaffold_origin_ids)
        }

        if (params.mode_upper == "CLIA") {
            // CLIA mode has no gamma/MLST/SRST2
            gamma_hv         = wtasmbld_report.map{ it -> create_empty_ch(it) }
            gamma_ar         = wtasmbld_report.map{ it -> create_empty_ch(it) }
            gamma_pf         = wtasmbld_report.map{ it -> create_empty_ch(it) }
            mlst             = wtasmbld_report.map{ it -> create_empty_ch(it) }
            fullgene_results = wtasmbld_report.map{ it -> create_empty_ch(it) }
        }

        if (params.mode_upper == "UPDATE_PHOENIX") {
            // Combining output based on id:meta.id to create pipeline stats file by sample -- is this verbose, ugly and annoying. yes, if anyone has a slicker way to do this we welcome the input. 
            pipeline_stats_ch = fastp_raw_qc.map{ meta, fastp_raw_qc           -> [[id:meta.id, project_id:meta.project_id],fastp_raw_qc]}\
                .join(fastp_total_qc.map{             meta, fastp_total_qc         -> [[id:meta.id, project_id:meta.project_id],fastp_total_qc]},         by: [0])\
                .join(fullgene_results.map{           meta, fullgene_results       -> [[id:meta.id, project_id:meta.project_id],fullgene_results]},       by: [0])\
                .join(trimd_report.map{               meta, report                 -> [[id:meta.id, project_id:meta.project_id],report]},                 by: [0])\
                .join(trimd_krona_html.map{           meta, trimd_krona_html       -> [[id:meta.id, project_id:meta.project_id],trimd_krona_html]},       by: [0])\
                .join(trimd_k2_bh_summary.map{        meta, trimd_k2_bh_summary    -> [[id:meta.id, project_id:meta.project_id],trimd_k2_bh_summary]},    by: [0])\
                .join(renamed_fastas.map{             meta, renamed_fastas         -> [[id:meta.id, project_id:meta.project_id],renamed_fastas]},         by: [0])\
                .join(filtered_fastas.map{            meta, filtered_fastas        -> [[id:meta.id, project_id:meta.project_id],filtered_fastas]},        by: [0])\
                .join(mlst.map{                       meta, mlst                   -> [[id:meta.id, project_id:meta.project_id],mlst]},                   by: [0])\
                .join(gamma_hv.map{                   meta, gamma_hv               -> [[id:meta.id, project_id:meta.project_id],gamma_hv]},               by: [0])\
                .join(gamma_ar.map{                   meta, gamma_ar               -> [[id:meta.id, project_id:meta.project_id],gamma_ar]},               by: [0])\
                .join(gamma_pf.map{                   meta, gamma_pf               -> [[id:meta.id, project_id:meta.project_id],gamma_pf]},               by: [0])\
                .join(quast_report.map{               meta, quast_report           -> [[id:meta.id, project_id:meta.project_id],quast_report]},           by: [0])\
                .join(busco.map{                      meta, busco                  -> [[id:meta.id, project_id:meta.project_id],busco]},                  by: [0])\
                .join(asmbld_report.map{              meta, asmbld_report          -> [[id:meta.id, project_id:meta.project_id],asmbld_report]},          by: [0])\
                .join(asmbld_krona_html.map{          meta, asmbld_krona_html      -> [[id:meta.id, project_id:meta.project_id],asmbld_krona_html]},      by: [0])\
                .join(asmbld_k2_bh_summary.map{       meta, asmbld_k2_bh_summary   -> [[id:meta.id, project_id:meta.project_id],asmbld_k2_bh_summary]},   by: [0])\
                .join(wtasmbld_krona_html.map{        meta, wtasmbld_krona_html    -> [[id:meta.id, project_id:meta.project_id],wtasmbld_krona_html]},    by: [0])\
                .join(wtasmbld_report.map{            meta, wtasmbld_report        -> [[id:meta.id, project_id:meta.project_id],wtasmbld_report]},        by: [0])\
                .join(wtasmbld_k2_bh_summary.map{     meta, wtasmbld_k2_bh_summary -> [[id:meta.id, project_id:meta.project_id],wtasmbld_k2_bh_summary]}, by: [0])\
                .join(taxa_id.map{                    meta, taxa_id                -> [[id:meta.id, project_id:meta.project_id],taxa_id]},                by: [0])\
                .join(format_ani.map{                 meta, format_ani             -> [[id:meta.id, project_id:meta.project_id],format_ani]},             by: [0])\
                .join(assembly_ratio.map{             meta, assembly_ratio         -> [[id:meta.id, project_id:meta.project_id],assembly_ratio]},         by: [0])\
                .join(amr_point_mutations.map{        meta, amr_point_mutations    -> [[id:meta.id, project_id:meta.project_id],amr_point_mutations]},    by: [0])\
                .join(gc_content.map{                 meta, gc_content             -> [[id:meta.id, project_id:meta.project_id],gc_content]},             by: [0])\
                .join(wtasmbld_report_with_rt.map{    meta, report, rt             -> [[id:meta.id, project_id:meta.project_id], rt.base] },              by: [0])
        } else {
            // Combining output based on id:meta.id to create pipeline stats file by sample -- is this verbose, ugly and annoying. yes, if anyone has a slicker way to do this we welcome the input. 
            pipeline_stats_ch = fastp_raw_qc.map{ meta, fastp_raw_qc           -> [[id:meta.id],fastp_raw_qc]}\
                .join(fastp_total_qc.map{             meta, fastp_total_qc         -> [[id:meta.id],fastp_total_qc]},         by: [0])\
                .join(fullgene_results.map{           meta, fullgene_results       -> [[id:meta.id],fullgene_results]},       by: [0])\
                .join(trimd_report.map{               meta, report                 -> [[id:meta.id],report]},                 by: [0])\
                .join(trimd_krona_html.map{           meta, trimd_krona_html       -> [[id:meta.id],trimd_krona_html]},       by: [0])\
                .join(trimd_k2_bh_summary.map{        meta, trimd_k2_bh_summary    -> [[id:meta.id],trimd_k2_bh_summary]},    by: [0])\
                .join(renamed_fastas.map{             meta, renamed_fastas         -> [[id:meta.id],renamed_fastas]},         by: [0])\
                .join(filtered_fastas.map{            meta, filtered_fastas        -> [[id:meta.id],filtered_fastas]},        by: [0])\
                .join(mlst.map{                       meta, mlst                   -> [[id:meta.id],mlst]},                   by: [0])\
                .join(gamma_hv.map{                   meta, gamma_hv               -> [[id:meta.id],gamma_hv]},               by: [0])\
                .join(gamma_ar.map{                   meta, gamma_ar               -> [[id:meta.id],gamma_ar]},               by: [0])\
                .join(gamma_pf.map{                   meta, gamma_pf               -> [[id:meta.id],gamma_pf]},               by: [0])\
                .join(quast_report.map{               meta, quast_report           -> [[id:meta.id],quast_report]},           by: [0])\
                .join(busco.map{                      meta, busco                  -> [[id:meta.id],busco]},                  by: [0])\
                .join(asmbld_report.map{              meta, asmbld_report          -> [[id:meta.id],asmbld_report]},          by: [0])\
                .join(asmbld_krona_html.map{          meta, asmbld_krona_html      -> [[id:meta.id],asmbld_krona_html]},      by: [0])\
                .join(asmbld_k2_bh_summary.map{       meta, asmbld_k2_bh_summary   -> [[id:meta.id],asmbld_k2_bh_summary]},   by: [0])\
                .join(wtasmbld_krona_html.map{        meta, wtasmbld_krona_html    -> [[id:meta.id],wtasmbld_krona_html]},    by: [0])\
                .join(wtasmbld_report.map{            meta, wtasmbld_report        -> [[id:meta.id],wtasmbld_report]},        by: [0])\
                .join(wtasmbld_k2_bh_summary.map{     meta, wtasmbld_k2_bh_summary -> [[id:meta.id],wtasmbld_k2_bh_summary]}, by: [0])\
                .join(taxa_id.map{                    meta, taxa_id                -> [[id:meta.id],taxa_id]},                by: [0])\
                .join(format_ani.map{                 meta, format_ani             -> [[id:meta.id],format_ani]},             by: [0])\
                .join(assembly_ratio.map{             meta, assembly_ratio         -> [[id:meta.id],assembly_ratio]},         by: [0])\
                .join(amr_point_mutations.map{        meta, amr_point_mutations    -> [[id:meta.id],amr_point_mutations]},    by: [0])\
                .join(gc_content.map{                 meta, gc_content             -> [[id:meta.id],gc_content]},             by: [0])\
                .join(wtasmbld_report_with_rt.map{    meta, report, rt             -> [[id:meta.id], rt.base] },              by: [0])
        }

        GENERATE_PIPELINE_STATS (
            pipeline_stats_ch, params.coverage
        )
        ch_versions = ch_versions.mix(GENERATE_PIPELINE_STATS.out.versions)

        pipeline_stats = GENERATE_PIPELINE_STATS.out.pipeline_stats


    emit:
        pipeline_stats  = pipeline_stats
        versions        = ch_versions // channel: [ versions.yml ]
}