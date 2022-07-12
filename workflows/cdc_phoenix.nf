/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPhoenix.initialise(params, log)


// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ] //removed , params.fasta to stop issue w/connecting to aws and igenomes not used
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters

//input on command line
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet/list not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES
========================================================================================
*/

include { ASSET_CHECK                    } from '../modules/local/asset_check'
include { RENAME_FASTA_HEADERS           } from '../modules/local/rename_fasta_headers'
include { BUSCO                          } from '../modules/local/busco'
include { GAMMA_S as GAMMA_PF            } from '../modules/local/gammas'
include { GAMMA as GAMMA_AR              } from '../modules/local/gamma'
include { GAMMA as GAMMA_HV              } from '../modules/local/gamma'
include { FASTP as FASTP_SINGLES         } from '../modules/local/fastp_singles'
include { BBMAP_REFORMAT                 } from '../modules/local/contig_less500'
include { QUAST                          } from '../modules/local/quast'
include { FASTANI                        } from '../modules/local/fastani'
include { DETERMINE_TOP_TAXA             } from '../modules/local/determine_top_taxa'
include { FORMAT_ANI                     } from '../modules/local/format_ANI_best_hit'
include { GATHERING_READ_QC_STATS        } from '../modules/local/fastp_minimizer'
include { DETERMINE_TAXA_ID              } from '../modules/local/tax_classifier'
include { CALCULATE_ASSEMBLY_RATIO       } from '../modules/local/assembly_ratio'
include { CREATE_SUMMARY_LINE            } from '../modules/local/phoenix_summary_line'
include { FETCH_FAILED_SUMMARIES         } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES           } from '../modules/local/phoenix_summary'
include { GENERATE_PIPELINE_STATS        } from '../modules/local/generate_pipeline_stats'

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK                    } from '../subworkflows/local/input_check'
include { SPADES_WF                      } from '../subworkflows/local/spades_failure'
include { GENERATE_PIPELINE_STATS_WF     } from '../subworkflows/local/generate_pipeline_stats'
include { KRAKEN2_WF as KRAKEN2_TRIMD    } from '../subworkflows/local/kraken2krona'
include { KRAKEN2_WF as KRAKEN2_ASMBLD   } from '../subworkflows/local/kraken2krona'
include { KRAKEN2_WF as KRAKEN2_WTASMBLD } from '../subworkflows/local/kraken2krona'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BBMAP_BBDUK                                             } from '../modules/nf-core/modules/bbmap/bbduk/main'
include { FASTP as FASTP_TRIMD                                    } from '../modules/nf-core/modules/fastp/main'
include { FASTQC as FASTQCTRIMD                                   } from '../modules/nf-core/modules/fastqc/main'
include { SRST2_SRST2 as SRST2_TRIMD_AR                           } from '../modules/nf-core/modules/srst2/srst2/main'
include { MLST                                                    } from '../modules/nf-core/modules/mlst/main'
include { MASH_DIST                                               } from '../modules/nf-core/modules/mash/dist/main'
include { MULTIQC                                                 } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                             } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { AMRFINDERPLUS_UPDATE                                    } from '../modules/nf-core/modules/amrfinderplus/update/main'
include { AMRFINDERPLUS_RUN                                       } from '../modules/nf-core/modules/amrfinderplus/run/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PHOENIX_EXQC {

    ch_versions     = Channel.empty() // Used to collect the software versions

    //
    // SUBWORKFLOW: Read in samplesheet/list, validate and stage input files
    //

    // Call in reads
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //unzip any zipped databases
    ASSET_CHECK (
        params.path2db
    )
    
    // Remove PhiX reads
    BBMAP_BBDUK (
        INPUT_CHECK.out.reads, params.bbdukdb
    )
    ch_versions = ch_versions.mix(BBMAP_BBDUK.out.versions)

    // Trim and remove low quality reads
    FASTP_TRIMD (
        BBMAP_BBDUK.out.reads, true, false
    )
    ch_versions = ch_versions.mix(FASTP_TRIMD.out.versions)

    // Rerun on unpaired reads to get stats, nothing removed
    FASTP_SINGLES (
        FASTP_TRIMD.out.reads_fail
    )
    ch_versions = ch_versions.mix(FASTP_SINGLES.out.versions)

    // Combining fastp json outputs based on meta.id
    fastp_json_ch = FASTP_TRIMD.out.json.join(FASTP_SINGLES.out.json, by: [0])

    // Script gathers data from jsons for pipeline stats file
    GATHERING_READ_QC_STATS(
        fastp_json_ch
    )

    // Running Fastqc on trimmed reads
    FASTQCTRIMD (
        FASTP_TRIMD.out.reads
    )
    ch_versions = ch_versions.mix(FASTQCTRIMD.out.versions.first())

    // Idenitifying AR genes in trimmed reads
    SRST2_TRIMD_AR (
        FASTP_TRIMD.out.reads.map{ meta, reads -> [ [id:meta.id, single_end:meta.single_end, db:'gene'], reads, params.ardb]}
    )
    ch_versions = ch_versions.mix(SRST2_TRIMD_AR.out.versions)

    // Checking for Contamination in trimmed reads, creating krona plots and best hit files
    KRAKEN2_TRIMD (
        FASTP_TRIMD.out.reads, "trimd", GATHERING_READ_QC_STATS.out.fastp_total_qc, []
    )
    ch_versions = ch_versions.mix(KRAKEN2_TRIMD.out.versions)

    SPADES_WF (
        FASTP_SINGLES.out.reads, FASTP_TRIMD.out.reads, \
        GATHERING_READ_QC_STATS.out.fastp_total_qc, \
        GATHERING_READ_QC_STATS.out.fastp_raw_qc, \
        SRST2_TRIMD_AR.out.fullgene_results, \
        KRAKEN2_TRIMD.out.report, KRAKEN2_TRIMD.out.krona_html, \
        KRAKEN2_TRIMD.out.k2_bh_summary, \
        true
    )
    ch_versions = ch_versions.mix(SPADES_WF.out.versions)

    // Rename scaffold headers
    RENAME_FASTA_HEADERS (
        SPADES_WF.out.spades_ch
    )
    ch_versions = ch_versions.mix(RENAME_FASTA_HEADERS.out.versions)

    // Removing scaffolds <500bp
    BBMAP_REFORMAT (
        RENAME_FASTA_HEADERS.out.renamed_scaffolds
    )
    ch_versions = ch_versions.mix(BBMAP_REFORMAT.out.versions)

    // Fetch AMRFinder Database
    AMRFINDERPLUS_UPDATE ( )
    ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
    
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

    // Checking for Contamination in assembly creating krona plots and best hit files
    KRAKEN2_ASMBLD (
        BBMAP_REFORMAT.out.reads,"asmbld", [], QUAST.out.report_tsv
    )
    ch_versions = ch_versions.mix(KRAKEN2_ASMBLD.out.versions)

    // Creating krona plots and best hit files for weighted assembly
    KRAKEN2_WTASMBLD (
        BBMAP_REFORMAT.out.reads,"wtasmbld", [], QUAST.out.report_tsv
    )
    ch_versions = ch_versions.mix(KRAKEN2_WTASMBLD.out.versions)

    // Running Mash distance to get top 20 matches for fastANI to speed things up
    MASH_DIST (
        BBMAP_REFORMAT.out.reads, params.mash_sketch
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
    best_hit_ch = KRAKEN2_WTASMBLD.out.report.map{meta, kraken_weighted_report -> [[id:meta.id], kraken_weighted_report]}\
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

    GENERATE_PIPELINE_STATS_WF (
        FASTP_TRIMD.out.reads, \
        GATHERING_READ_QC_STATS.out.fastp_raw_qc, \
        GATHERING_READ_QC_STATS.out.fastp_total_qc, \
        SRST2_TRIMD_AR.out.fullgene_results, \
        KRAKEN2_TRIMD.out.report, \
        KRAKEN2_TRIMD.out.krona_html, \
        KRAKEN2_TRIMD.out.k2_bh_summary, \
        RENAME_FASTA_HEADERS.out.renamed_scaffolds, \
        BBMAP_REFORMAT.out.reads, \
        MLST.out.tsv, \
        GAMMA_HV.out.gamma, \
        GAMMA_AR.out.gamma, \
        GAMMA_PF.out.gamma, \
        QUAST.out.report_tsv, \
        BUSCO.out.short_summaries_specific_txt, \
        KRAKEN2_ASMBLD.out.report, \
        KRAKEN2_ASMBLD.out.krona_html, \
        KRAKEN2_ASMBLD.out.k2_bh_summary, \
        KRAKEN2_WTASMBLD.out.krona_html, \
        KRAKEN2_WTASMBLD.out.report, \
        KRAKEN2_WTASMBLD.out.k2_bh_summary, \
        DETERMINE_TAXA_ID.out.taxonomy, \
        FORMAT_ANI.out.ani_best_hit, \
        CALCULATE_ASSEMBLY_RATIO.out.ratio,\
        true
    )

    // Combining output based on meta.id to create summary by sample -- is this verbose, ugly and annoying? yes, if anyone has a slicker way to do this we welcome the input. 
    line_summary_ch = GATHERING_READ_QC_STATS.out.fastp_total_qc.map{meta, fastp_total_qc -> [[id:meta.id], fastp_total_qc]}\
    .join(MLST.out.tsv.map{                                          meta, tsv            -> [[id:meta.id], tsv]},            by: [0])\
    .join(GAMMA_HV.out.gamma.map{                                    meta, gamma          -> [[id:meta.id], gamma]},          by: [0])\
    .join(GAMMA_AR.out.gamma.map{                                    meta, gamma          -> [[id:meta.id], gamma]},          by: [0])\
    .join(QUAST.out.report_tsv.map{                                  meta, report_tsv     -> [[id:meta.id], report_tsv]},     by: [0])\
    .join(CALCULATE_ASSEMBLY_RATIO.out.ratio.map{                    meta, ratio          -> [[id:meta.id], ratio]},          by: [0])\
    .join(GENERATE_PIPELINE_STATS_WF.out.pipeline_stats.map{         meta, pipeline_stats -> [[id:meta.id], pipeline_stats]}, by: [0])\
    .join(DETERMINE_TAXA_ID.out.taxonomy.map{                        meta, taxonomy       -> [[id:meta.id], taxonomy]},       by: [0])\
    .join(KRAKEN2_TRIMD.out.k2_bh_summary.map{                       meta, k2_bh_summary  -> [[id:meta.id], k2_bh_summary]},  by: [0])

    // Generate summary per sample
    CREATE_SUMMARY_LINE(
        line_summary_ch
    )
    ch_versions = ch_versions.mix(CREATE_SUMMARY_LINE.out.versions)

    // Collect all the summary files prior to fetch step to force the fetch process to wait to wait
    failed_summaries_ch         = SPADES_WF.out.line_summary.collect()
    summaries_ch                = CREATE_SUMMARY_LINE.out.line_summary.collect()

    // combine all line summaries into one channel
    FETCH_FAILED_SUMMARIES (
        params.outdir, failed_summaries_ch, summaries_ch
    )

    spades_failure_summaries_ch = FETCH_FAILED_SUMMARIES.out.spades_failure_summary_line
    all_summaries_ch = spades_failure_summaries_ch.combine(failed_summaries_ch).combine(summaries_ch)

    // Combining sample summaries into final report
    GATHER_SUMMARY_LINES (
        all_summaries_ch
    )
    ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

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

/*
========================================================================================
    THE END
========================================================================================
*/
