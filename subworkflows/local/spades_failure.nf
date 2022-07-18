//
// Subworkflow: Running SPAdes and checking if spades failed to create scaffolds
//

include { SPADES                               } from '../../modules/local/spades'
include { CREATE_SUMMARY_LINE_FAILURE          } from '../../modules/local/phoenix_summary_line_failure'
include { GENERATE_PIPELINE_STATS_FAILURE      } from '../../modules/local/generate_pipeline_stats_failure'
include { GENERATE_PIPELINE_STATS_FAILURE_EXQC } from '../../modules/local/generate_pipeline_stats_failure_exqc'
include { DETERMINE_TAXA_ID_FAILURE            } from '../../modules/local/taxa_classifier_failure'

workflow SPADES_WF {
    take:
        single_reads      // channel: tuple val(meta), path(reads), path(single_reads): FASTP_SINGLES.out.reads
        paired_reads      // channel: tuple val(meta), path(reads), path(paired_reads):FASTP_TRIMD.out.reads.map
        fastp_total_qc    // channel: tuple (meta) path(fastp_total_qc): GATHERING_READ_QC_STATS.out.fastp_total_qc
        fastp_raw_qc      // channel: tuple (meta) path(fastp_raw_qc): GATHERING_READ_QC_STATS.out.fastp_raw_qc
        fullgene_results  // channel: tuple (meta) path(fullgene_results): SRST2_TRIMD_AR.out.fullgene_results
        report            // channel: tuple (meta) path(report): KRAKEN2_TRIMD.out.report
        krona_html        // channel: tuple (meta) path(krona_html): KRAKEN2_TRIMD.out.krona_html
        k2_bh_summary     // channel: tuple (meta) path(k2_bh_summary): KRAKEN2_TRIMD.out.k2_bh_summary
        extended_qc

    main:
        ch_versions     = Channel.empty() // Used to collect the software versions

        // Combining paired end reads and unpaired reads that pass QC filters, both get passed to Spades
        passing_reads_ch = paired_reads.map{ meta, reads          -> [[id:meta.id],reads]}\
        .join(single_reads.map{              meta, reads          -> [[id:meta.id],reads]},          by: [0])\
        .join(k2_bh_summary.map{             meta, ksummary       -> [[id:meta.id],ksummary]},       by: [0])\
        .join(fastp_raw_qc.map{              meta, fastp_raw_qc   -> [[id:meta.id],fastp_raw_qc]},   by: [0])\
        .join(fastp_total_qc.map{            meta, fastp_total_qc -> [[id:meta.id],fastp_total_qc]}, by: [0])\
        .join(report.map{                    meta, report         -> [[id:meta.id],report]},         by: [0])\
        .join(krona_html.map{                meta, krona_html     -> [[id:meta.id],krona_html]},     by: [0])

        // Assemblying into scaffolds by passing filtered paired in reads and unpaired reads
        SPADES (
            passing_reads_ch
        )
        ch_versions = ch_versions.mix(SPADES.out.versions)

        if (extended_qc == true) { // Run this if there extra stats are requested via --extended_qc

            // Combining weighted kraken report with the FastANI hit based on meta.id
            best_hit_ch = k2_bh_summary.map{                         meta, ksummary       -> [[id:meta.id], ksummary]}\
            .join(SPADES.out.spades_outcome.splitCsv(strip:true).map{meta, spades_outcome -> [[id:meta.id], spades_outcome]})

            // Getting ID from either FastANI or if fails, from Kraken2
            DETERMINE_TAXA_ID_FAILURE (
                    best_hit_ch, params.taxa
            )

            pipeline_stats_ch = paired_reads.map{            meta, reads            -> [[id:meta.id],reads]}\
            .join(fastp_raw_qc.map{                          meta, fastp_raw_qc     -> [[id:meta.id],fastp_raw_qc]},     by: [0])\
            .join(fastp_total_qc.map{                        meta, fastp_total_qc   -> [[id:meta.id],fastp_total_qc]},   by: [0])\
            .join(fullgene_results.map{                      meta, fullgene_results -> [[id:meta.id],fullgene_results]}, by: [0])\
            .join(report.map{                                meta, report           -> [[id:meta.id],report]},           by: [0])\
            .join(krona_html.map{                            meta, html             -> [[id:meta.id],html]},             by: [0])\
            .join(k2_bh_summary.map{                         meta, ksummary         -> [[id:meta.id],ksummary]},         by: [0])\
            .join(DETERMINE_TAXA_ID_FAILURE.out.taxonomy.map{meta, taxonomy         -> [[id:meta.id],taxonomy]},         by: [0])
            
            // Adding the outcome of spades (scaffolds created or not) to the channel
            pipeline_stats_ch = pipeline_stats_ch.join(SPADES.out.spades_outcome.splitCsv(strip:true).map{meta, spades_outcome -> [[id:meta.id], spades_outcome]})
            
            // Generate pipeline stats for case when spades fails to create scaffolds
            GENERATE_PIPELINE_STATS_FAILURE_EXQC (
                pipeline_stats_ch
            )

            // Adding in trimmed reads info into channel
            line_summary_ch = GENERATE_PIPELINE_STATS_FAILURE_EXQC.out.pipeline_stats.map{ meta, pipeline_stats  -> [[id:meta.id],pipeline_stats]}\
            .join(fastp_total_qc.map{                                                      meta, fastp_total_qc  -> [[id:meta.id],fastp_total_qc]}, by: [0])\
            .join(k2_bh_summary.map{                                                       meta, ksummary        -> [[id:meta.id],ksummary]},       by: [0])\
            .join(DETERMINE_TAXA_ID_FAILURE.out.taxonomy.map{meta, taxonomy -> [[id:meta.id],taxonomy]},         by: [0])\

            // Adding the outcome of spades (scaffolds created or not) to the channel 
            line_summary_ch = line_summary_ch.join(SPADES.out.spades_outcome.splitCsv(strip:true).map{meta, spades_outcome -> [[id:meta.id], spades_outcome]})

        } else {

            // Combining weighted kraken report with the FastANI hit based on meta.id
            best_hit_ch = k2_bh_summary.map{                         meta, ksummary       -> [[id:meta.id], ksummary]}\
            .join(SPADES.out.spades_outcome.splitCsv(strip:true).map{meta, spades_outcome -> [[id:meta.id], spades_outcome]})

            // Getting ID from either FastANI or if fails, from Kraken2
            DETERMINE_TAXA_ID_FAILURE (
                best_hit_ch, params.taxa
            )

            pipeline_stats_ch = paired_reads.map{            meta, reads           -> [[id:meta.id],reads]}\
            .join(fastp_raw_qc.map{                          meta, fastp_raw_qc    -> [[id:meta.id],fastp_raw_qc]},     by: [0])\
            .join(fastp_total_qc.map{                        meta, fastp_total_qc  -> [[id:meta.id],fastp_total_qc]},   by: [0])\
            .join(report.map{                                meta, report          -> [[id:meta.id],report]},           by: [0])\
            .join(krona_html.map{                            meta, html            -> [[id:meta.id],html]},             by: [0])\
            .join(k2_bh_summary.map{                         meta, ksummary        -> [[id:meta.id],ksummary]},         by: [0])\
            .join(DETERMINE_TAXA_ID_FAILURE.out.taxonomy.map{meta, taxonomy        -> [[id:meta.id],taxonomy]},         by: [0])\

            // Adding the outcome of spades (scaffolds created or not) to the channel
            pipeline_stats_ch = pipeline_stats_ch.join(SPADES.out.spades_outcome.splitCsv(strip:true).map{meta, spades_outcome -> [[id:meta.id], spades_outcome]})

            // Generate pipeline stats for case when spades fails to create scaffolds
            GENERATE_PIPELINE_STATS_FAILURE (
                pipeline_stats_ch
            )

            // Adding in trimmed reads info into channel
            line_summary_ch = GENERATE_PIPELINE_STATS_FAILURE.out.pipeline_stats.map{ meta, pipeline_stats  -> [[id:meta.id],pipeline_stats]}\
            .join(fastp_total_qc.map{                                                 meta, fastp_total_qc  -> [[id:meta.id],fastp_total_qc]}, by: [0])\
            .join(k2_bh_summary.map{                                                  meta, ksummary        -> [[id:meta.id],ksummary]},       by: [0])
            .join(DETERMINE_TAXA_ID_FAILURE.out.taxonomy.map{meta, taxonomy -> [[id:meta.id],taxonomy]},         by: [0])\

            // Adding the outcome of spades (scaffolds created or not) to the channel
            line_summary_ch = line_summary_ch.join(SPADES.out.spades_outcome.splitCsv(strip:true).map{meta, spades_outcome -> [[id:meta.id], spades_outcome]})
        }

        // Create one line summary for case when spades fails to create scaffolds
        CREATE_SUMMARY_LINE_FAILURE (
            line_summary_ch
        )

        // Defining out channel
        spades_ch = SPADES.out.scaffolds.map{meta, scaffolds -> [ [id:meta.id, single_end:true], scaffolds]}

    emit:
        spades_ch                   = spades_ch
        spades_outcome              = SPADES.out.spades_outcome
        line_summary                = CREATE_SUMMARY_LINE_FAILURE.out.line_summary
        versions                    = ch_versions // channel: [ versions.yml ]
}