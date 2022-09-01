include { GAMMA } from '../modules/nf-core/modules/gamma/main'
include { FASTP } from '../modules/nf-core/modules/fastp/main'
include { SPADES } from '../modules/nf-core/modules/spades/main'
include { QUAST } from '../modules/nf-core/modules/quast/main'
include { FASTANI } from '../modules/nf-core/modules/fastani/main'
include { MASH_DIST } from '../modules/nf-core/modules/mash/dist/main'
include { MLST } from '../modules/nf-core/modules/mlst/main'
include { KRAKEN2 } from '../modules/nf-core/modules/kraken2/main'
include { PROKKA } from '../modules/nf-core/modules/prokka/main'
include { BUSCO } from '../modules/nf-core/modules/busco/main'
include { KRONA } from '../modules/nf-core/modules/krona/main'

//local module
include { KRAKEN2_DB } from '../modules/local/kraken2db'

// database parameter checks

if(params.gamma_db){
    Channel
        .fromPath( "${params.gamma_db}" )
        .set { ch_gamma }
} else {
    ch_gamma = Channel.empty()
}


workflow SPADES_TO_END {
    
    // Assemblying into scaffolds by passing filtered paired in reads and unpaired reads
    SPADES_LOCAL (
        passing_reads_ch
    )
    ch_versions = ch_versions.mix(SPADES_LOCAL.out.versions)
    spades_ch = SPADES_LOCAL.out.scaffolds.map{meta, scaffolds -> [ [id:meta.id, single_end:true], scaffolds]}

    // Rename scaffold headers
    RENAME_FASTA_HEADERS (
        spades_ch
    )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADERS.out.versions)

    // Removing scaffolds <500bp
    BBMAP_REFORMAT (
        RENAME_FASTA_HEADERS.out.renamed_scaffolds
    )
    ch_versions = ch_versions.mix(BBMAP_REFORMAT.out.versions)

    // Getting MLST scheme for taxa
    MLST (
        BBMAP_REFORMAT.out.reads
    )
    ch_versions = ch_versions.mix(MLST.out.versions)

    // Running gamma to identify hypervirulence genes in scaffolds
    GAMMA_HV (
        BBMAP_REFORMAT.out.reads, params.hvgamdb
    )
    ch_versions = ch_versions.mix(GAMMA_HV.out.versions)

    // Running gamma to identify AR genes in scaffolds
    GAMMA_AR (
        BBMAP_REFORMAT.out.reads, params.ardb
    )
    ch_versions = ch_versions.mix(GAMMA_AR.out.versions)

    GAMMA_PF (
        BBMAP_REFORMAT.out.reads, params.gamdbpf
    )
    ch_versions = ch_versions.mix(GAMMA_PF.out.versions)

    // Getting Assembly Stats
    QUAST (
        BBMAP_REFORMAT.out.reads
    )
    ch_versions = ch_versions.mix(QUAST.out.versions)

    // Checking single copy genes for assembly completeness
    BUSCO (
        BBMAP_REFORMAT.out.reads, 'auto', [], []
    )
    ch_versions = ch_versions.mix(BUSCO.out.versions)

    // Getting species ID as back up for FastANI and checking contamination isn't in assembly
    KRAKEN2_ASMBLD (
        BBMAP_REFORMAT.out.reads, params.path2db, "asmbld", true, true
    )
    ch_versions = ch_versions.mix(KRAKEN2_ASMBLD.out.versions)

    // Create mpa file
    KREPORT2MPA_ASMBLD (
        KRAKEN2_ASMBLD.out.report
    )
    ch_versions = ch_versions.mix(KREPORT2MPA_ASMBLD.out.versions)

    // Converting kraken report to krona file to have hierarchical output in krona plot
    KREPORT2KRONA_ASMBLD (
        KRAKEN2_ASMBLD.out.report, "asmbld"
    )
    ch_versions = ch_versions.mix(KREPORT2KRONA_ASMBLD.out.versions)

    // Create krona plot from kraken report 
    KRONA_KTIMPORTTEXT_ASMBLD (
        KREPORT2KRONA_ASMBLD.out.krona, "asmbld"
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT_ASMBLD.out.versions)

    // Combining kraken report with quast report based on meta.id
    kraken_bh_asmbld_ch = KRAKEN2_ASMBLD.out.report.map{meta, report     -> [[id:meta.id], report]}\
    .join(QUAST.out.report_tsv.map{                     meta, report_tsv -> [[id:meta.id], report_tsv]}, by: [0])

    // Getting Kraken best hit for assembled data
    KRAKEN2_BH_ASMBLD (
        kraken_bh_asmbld_ch, "asmbld"
    )

    // Getting species ID as back up for FastANI and checking contamination isn't in assembly
    KRAKEN2_ASMBLD_WEIGHTED (
        BBMAP_REFORMAT.out.reads, params.path2db, "wtasmbld", true, true
    )
    ch_versions = ch_versions.mix(KRAKEN2_ASMBLD_WEIGHTED.out.versions)

    // Create weighted kraken report based on scaffold length
    KRAKENTOOLS_MAKEKREPORT (
        KRAKEN2_ASMBLD_WEIGHTED.out.classified_reads_assignment, params.ktaxmap
    )
    ch_versions = ch_versions.mix(KRAKENTOOLS_MAKEKREPORT.out.versions)

    // Converting kraken report to krona file to have hierarchical output in krona plot
    KREPORT2KRONA_WTASMBLD (
        KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report, "wtasmbld"
    )
    ch_versions = ch_versions.mix(KREPORT2KRONA_WTASMBLD.out.versions)

    // Combining kraken report with quast report based on meta.id
    kraken_bh_wtasmbld_ch = KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report.map{meta, kraken_weighted_report -> [[id:meta.id], kraken_weighted_report]}\
    .join(QUAST.out.report_tsv.map{                                                meta, report_tsv             -> [[id:meta.id], report_tsv]}, by: [0])

    // Getting Kraken best hit for assembled data
    KRAKEN2_BH_ASMBLD_WEIGHTED(
        kraken_bh_wtasmbld_ch, "wtasmbld"
    )

    KRONA_KTIMPORTTEXT_WTASMBLD (
        KREPORT2KRONA_WTASMBLD.out.krona, "wtasmbld"
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT_WTASMBLD.out.versions)

    // Running Mash distance to get top 20 matches for fastANI to speed things up
    MASH_DIST (
        BBMAP_REFORMAT.out.reads, ASSET_CHECK.out.mash_sketch
    )
    ch_versions = ch_versions.mix(MASH_DIST.out.versions)

    // Combining mash dist with filtered scaffolds based on meta.id
    top_taxa_ch = MASH_DIST.out.dist.map{ meta, dist  -> [[id:meta.id], dist]}\
    .join(BBMAP_REFORMAT.out.reads.map{   meta, reads -> [[id:meta.id], reads ]}, by: [0])

    // Generate file with list of paths of top taxa for fastANI
    DETERMINE_TOP_TAXA (
        top_taxa_ch, params.refseq_fasta_database
    )

    // Combining filtered scaffolds with the top taxa list based on meta.id
    top_taxa_list_ch = BBMAP_REFORMAT.out.reads.map{meta, reads         -> [[id:meta.id], reads]}\
    .join(DETERMINE_TOP_TAXA.out.top_taxa_list.map{ meta, top_taxa_list -> [[id:meta.id], top_taxa_list ]}, by: [0])

    // Getting species ID
    FASTANI (
        top_taxa_list_ch, params.refseq_fasta_database
    )
    ch_versions = ch_versions.mix(FASTANI.out.versions)

    // Reformat ANI headers
    FORMAT_ANI (
        FASTANI.out.ani
    )

    // Combining weighted kraken report with the FastANI hit based on meta.id
    best_hit_ch = KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report.map{meta, kraken_weighted_report -> [[id:meta.id], kraken_weighted_report]}\
    .join(FORMAT_ANI.out.ani_best_hit.map{                               meta, ani_best_hit           -> [[id:meta.id], ani_best_hit ]}, by: [0])

    // Getting ID from either FastANI or if fails, from Kraken2
    DETERMINE_TAXA_ID (
        best_hit_ch, params.taxa
    )

    // Combining determined taxa with the assembly stats based on meta.id
    assembly_ratios_ch = DETERMINE_TAXA_ID.out.taxonomy.map{meta, taxonomy   -> [[id:meta.id], taxonomy]}\
    .join(QUAST.out.report_tsv.map{                         meta, report_tsv -> [[id:meta.id], report_tsv]}, by: [0])

    // Calculating the assembly ratio
    CALCULATE_ASSEMBLY_RATIO (
        assembly_ratios_ch, params.ncbi_assembly_stats
    )

    // Combining output based on meta.id to create summary by sample -- is this verbose, ugly and annoying? yes, if anyone has a slicker way to do this we welcome the input. 
    line_summary_ch = GATHERING_READ_QC_STATS.out.fastp_total_qc.map{meta, fastp_total_qc -> [[id:meta.id], fastp_total_qc]}\
    .join(MLST.out.tsv.map{                                          meta, tsv            -> [[id:meta.id], tsv]},        by: [0])\
    .join(GAMMA_HV.out.gamma.map{                                    meta, gamma          -> [[id:meta.id], gamma]},      by: [0])\
    .join(GAMMA_AR.out.gamma.map{                                    meta, gamma          -> [[id:meta.id], gamma]},      by: [0])\
    .join(QUAST.out.report_tsv.map{                                  meta, report_tsv     -> [[id:meta.id], report_tsv]}, by: [0])\
    .join(CALCULATE_ASSEMBLY_RATIO.out.ratio.map{                    meta, ratio          -> [[id:meta.id], ratio]},      by: [0])

    // Generate summary per sample
    CREATE_SUMMARY_LINE(
        line_summary_ch
    )
    ch_versions = ch_versions.mix(CREATE_SUMMARY_LINE.out.versions)

    // combine all line summaries into one channel
    all_summaries_ch = CREATE_SUMMARY_LINE.out.line_summary.collect()

    // Combining sample summaries into final report
    GATHER_SUMMARY_LINES (
        all_summaries_ch
    )
    ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

    // Combining output based on id:meta.id to create pipeline stats file by sample -- is this verbose, ugly and annoying. yes, if anyone has a slicker way to do this we welcome the input. 
    pipeline_stats_ch = FASTP_TRIMD.out.reads.map{                    meta, reads                        -> [[id:meta.id],reads]}\
    .join(GATHERING_READ_QC_STATS.out.fastp_raw_qc.map{               meta, fastp_raw_qc                 -> [[id:meta.id],fastp_raw_qc]},                 by: [0])\
    .join(GATHERING_READ_QC_STATS.out.fastp_total_qc.map{             meta, fastp_total_qc               -> [[id:meta.id],fastp_total_qc]},               by: [0])\
    .join(SRST2_TRIMD_AR.out.fullgene_results.map{                    meta, fullgene_results             -> [[id:meta.id],fullgene_results]},             by: [0])\
    .join(KRAKEN2_TRIMD.out.report.map{                               meta, report                       -> [[id:meta.id],report]},                       by: [0])\
    .join(KRONA_KTIMPORTTEXT_TRIMD.out.html.map{                      meta, html                         -> [[id:meta.id],html]},                         by: [0])\
    .join(KRAKEN2_BH_TRIMD.out.ksummary.map{                          meta, ksummary                     -> [[id:meta.id],ksummary]},                     by: [0])\
    .join(RENAME_FASTA_HEADERS.out.renamed_scaffolds.map{             meta, renamed_scaffolds            -> [[id:meta.id],renamed_scaffolds]},            by: [0])\
    .join(BBMAP_REFORMAT.out.reads.map{                               meta, reads                        -> [[id:meta.id],reads]},                        by: [0])\
    .join(MLST.out.tsv.map{                                           meta, tsv                          -> [[id:meta.id],tsv]},                          by: [0])\
    .join(GAMMA_HV.out.gamma.map{                                     meta, gamma                        -> [[id:meta.id],gamma]},                        by: [0])\
    .join(GAMMA_AR.out.gamma.map{                                     meta, gamma                        -> [[id:meta.id],gamma]},                        by: [0])\
    .join(GAMMA_PF.out.gamma.map{                                     meta, gamma                        -> [[id:meta.id],gamma]},                        by: [0])\
    .join(QUAST.out.report_tsv.map{                                   meta, report_tsv                   -> [[id:meta.id],report_tsv]},                   by: [0])\
    .join(BUSCO.out.short_summaries_specific_txt.map{                 meta, short_summaries_specific_txt -> [[id:meta.id],short_summaries_specific_txt]}, by: [0])\
    .join(KRAKEN2_ASMBLD.out.report.map{                              meta, report                       -> [[id:meta.id],report]},                       by: [0])\
    .join(KRONA_KTIMPORTTEXT_ASMBLD.out.html.map{                     meta, html                         -> [[id:meta.id],html]},                         by: [0])\
    .join(KRAKEN2_BH_ASMBLD.out.ksummary.map{                         meta, ksummary                     -> [[id:meta.id],ksummary]},                     by: [0])\
    .join(KRONA_KTIMPORTTEXT_WTASMBLD.out.html.map{                   meta, html                         -> [[id:meta.id],html]},                         by: [0])\
    .join(KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report.map{     meta, kraken_weighted_report       -> [[id:meta.id],kraken_weighted_report]},       by: [0])\
    .join(KRAKEN2_BH_ASMBLD_WEIGHTED.out.ksummary.map{                meta, ksummary                     -> [[id:meta.id],ksummary]},                     by: [0])\
    .join(DETERMINE_TAXA_ID.out.taxonomy.map{                         meta, taxonomy                     -> [[id:meta.id],taxonomy]},                     by: [0])\
    .join(FORMAT_ANI.out.ani_best_hit.map{                            meta, ani_best_hit                 -> [[id:meta.id],ani_best_hit]},                 by: [0])\
    .join(CALCULATE_ASSEMBLY_RATIO.out.ratio.map{                     meta, ratio                        -> [[id:meta.id],ratio]},                        by: [0])

    GENERATE_PIPELINE_STATS (
        pipeline_stats_ch
    )

    // Collecting the software versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowPhoenix.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQCTRIMD.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}