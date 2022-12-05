/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowPhoenix.initialise(params, log)


// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db] //removed , params.fasta to stop issue w/connecting to aws and igenomes not used
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters

//input on command line
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet/list not specified!' }
if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified!' }

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

include { RENAME_FASTA_HEADERS           } from '../modules/local/rename_fasta_headers'
include { GAMMA_S as GAMMA_PF            } from '../modules/local/gammas'
include { GAMMA as GAMMA_AR              } from '../modules/local/gamma'
include { GAMMA as GAMMA_HV              } from '../modules/local/gamma'
include { MLST                           } from '../modules/local/mlst'
include { BBMAP_REFORMAT                 } from '../modules/local/contig_less500'
include { QUAST                          } from '../modules/local/quast'
include { MASH_DIST                      } from '../modules/local/mash_distance'
include { FASTANI                        } from '../modules/local/fastani'
include { DETERMINE_TOP_TAXA             } from '../modules/local/determine_top_taxa'
include { FORMAT_ANI                     } from '../modules/local/format_ANI_best_hit'
include { GATHERING_READ_QC_STATS        } from '../modules/local/fastp_minimizer'
include { DETERMINE_TAXA_ID              } from '../modules/local/tax_classifier'
include { PROKKA                         } from '../modules/local/prokka'
include { AMRFINDERPLUS_UPDATE           } from '../modules/local/update_amrfinder_db'
include { GET_TAXA_FOR_AMRFINDER         } from '../modules/local/get_taxa_for_amrfinder'
include { AMRFINDERPLUS_RUN              } from '../modules/local/run_amrfinder'
include { CALCULATE_ASSEMBLY_RATIO       } from '../modules/local/assembly_ratio'
include { CREATE_SUMMARY_LINE            } from '../modules/local/phoenix_summary_line'
include { GATHER_SUMMARY_LINES           } from '../modules/local/phoenix_summary'

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/

include { SCAFFOLDS_INPUT_CHECK          } from '../subworkflows/local/scaffolds_input_check'
include { GENERATE_PIPELINE_STATS_WF     } from '../subworkflows/local/generate_pipeline_stats'
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

include { MULTIQC                      } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []
def count = 0

workflow SCAFFOLD_EXTERNAL {
    
        ch_versions     = Channel.empty() // Used to collect the software versions
        
        //Create samplesheet
        SCAFFOLDS_INPUT_CHECK (
            params.scaffolds_samplesheet
        )
        ch_versions = ch_versions.mix(SCAFFOLDS_INPUT_CHECK.out.versions)
        
        //spades_ch = SCAFFOLDS_INPUT_CHECK.out.scaffolds.map{meta, scaffolds -> [ [id:meta.id, single_end:true], scaffolds]}
        
        // Rename scaffold headers
        RENAME_FASTA_HEADERS (
            SCAFFOLDS_INPUT_CHECK.out.scaffolds
        )
        ch_versions = ch_versions.mix(RENAME_FASTA_HEADERS.out.versions)

        // Removing scaffolds <500bp
        BBMAP_REFORMAT (
            RENAME_FASTA_HEADERS.out.renamed_scaffolds
        )
        ch_versions = ch_versions.mix(BBMAP_REFORMAT.out.versions)

        // Getting MLST scheme for taxa
        MLST (
            BBMAP_REFORMAT.out.filtered_scaffolds
        )
        ch_versions = ch_versions.mix(MLST.out.versions)

        // Running gamma to identify hypervirulence genes in scaffolds
        GAMMA_HV (
            BBMAP_REFORMAT.out.filtered_scaffolds, params.hvgamdb
        )
        ch_versions = ch_versions.mix(GAMMA_HV.out.versions)

        // Running gamma to identify AR genes in scaffolds
        GAMMA_AR (
            BBMAP_REFORMAT.out.filtered_scaffolds, params.ardb
        )
        ch_versions = ch_versions.mix(GAMMA_AR.out.versions)

        GAMMA_PF (
            BBMAP_REFORMAT.out.filtered_scaffolds, params.gamdbpf
        )
        ch_versions = ch_versions.mix(GAMMA_PF.out.versions)

        // Getting Assembly Stats
        QUAST (
            BBMAP_REFORMAT.out.filtered_scaffolds
        )
        ch_versions = ch_versions.mix(QUAST.out.versions)

        // Creating krona plots and best hit files for weighted assembly
        KRAKEN2_WTASMBLD (
            BBMAP_REFORMAT.out.filtered_scaffolds,"wtasmbld", [], QUAST.out.report_tsv
        )
        ch_versions = ch_versions.mix(KRAKEN2_WTASMBLD.out.versions)

        // Running Mash distance to get top 20 matches for fastANI to speed things up
        MASH_DIST (
            BBMAP_REFORMAT.out.filtered_scaffolds, ASSET_CHECK.out.mash_sketch
        )
        ch_versions = ch_versions.mix(MASH_DIST.out.versions)

        // Combining mash dist with filtered scaffolds based on meta.id
        top_taxa_ch = MASH_DIST.out.dist.map{ meta, dist  -> [[id:meta.id], dist]}\
        .join(BBMAP_REFORMAT.out.filtered_scaffolds.map{   meta, reads -> [[id:meta.id], reads ]}, by: [0])

        // Generate file with list of paths of top taxa for fastANI
        DETERMINE_TOP_TAXA (
            top_taxa_ch
        )

        // Combining filtered scaffolds with the top taxa list based on meta.id
        top_taxa_list_ch = BBMAP_REFORMAT.out.filtered_scaffolds.map{ meta, reads           -> [[id:meta.id], reads]}\
        .join(DETERMINE_TOP_TAXA.out.top_taxa_list.map{  meta, top_taxa_list   -> [[id:meta.id], top_taxa_list ]}, by: [0])\
        .join(DETERMINE_TOP_TAXA.out.reference_files.map{meta, reference_files -> [[id:meta.id], reference_files ]}, by: [0])

        // Getting species ID
        FASTANI (
            top_taxa_list_ch
        )
        ch_versions = ch_versions.mix(FASTANI.out.versions)

        // Reformat ANI headers
        FORMAT_ANI (
            FASTANI.out.ani
        )

        // Combining weighted kraken report with the FastANI hit based on meta.id
        best_hit_ch = KRAKEN2_WTASMBLD.out.report.map{meta, kraken_weighted_report -> [[id:meta.id], kraken_weighted_report]}\
        .join(FORMAT_ANI.out.ani_best_hit.map{        meta, ani_best_hit           -> [[id:meta.id], ani_best_hit ]},  by: [0])\
        .join(KRAKEN2_TRIMD.out.k2_bh_summary.map{    meta, k2_bh_summary          -> [[id:meta.id], k2_bh_summary ]}, by: [0])

        // Getting ID from either FastANI or if fails, from Kraken2
        DETERMINE_TAXA_ID (
            best_hit_ch, params.taxa
        )
        ch_versions = ch_versions.mix(DETERMINE_TAXA_ID.out.versions)

        // get gff and protein files for amrfinder+
        PROKKA (
            BBMAP_REFORMAT.out.filtered_scaffolds, [], []
        )
        ch_versions = ch_versions.mix(PROKKA.out.versions)

        // Fetch AMRFinder Database
        AMRFINDERPLUS_UPDATE( )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)

        // Create file that has the organism name to pass to AMRFinder
        GET_TAXA_FOR_AMRFINDER (
            DETERMINE_TAXA_ID.out.taxonomy
        )

        // Combining taxa and scaffolds to run amrfinder and get the point mutations. 
        amr_channel = BBMAP_REFORMAT.out.filtered_scaffolds.map{                               meta, reads          -> [[id:meta.id], reads]}\
        .join(GET_TAXA_FOR_AMRFINDER.out.amrfinder_taxa.splitCsv(strip:true).map{ meta, amrfinder_taxa -> [[id:meta.id], amrfinder_taxa ]}, by: [0])\
        .join(PROKKA.out.faa.map{                                                 meta, faa            -> [[id:meta.id], faa ]},            by: [0])\
        .join(PROKKA.out.gff.map{                                                 meta, gff            -> [[id:meta.id], gff ]},            by: [0])

        // Run AMRFinder
        AMRFINDERPLUS_RUN (
            amr_channel, AMRFINDERPLUS_UPDATE.out.db
        )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

        // Combining determined taxa with the assembly stats based on meta.id
        assembly_ratios_ch = DETERMINE_TAXA_ID.out.taxonomy.map{meta, taxonomy   -> [[id:meta.id], taxonomy]}\
        .join(QUAST.out.report_tsv.map{                         meta, report_tsv -> [[id:meta.id], report_tsv]}, by: [0])

        // Calculating the assembly ratio and gather GC% stats
        CALCULATE_ASSEMBLY_RATIO (
            assembly_ratios_ch, params.ncbi_assembly_stats
        )
        ch_versions = ch_versions.mix(CALCULATE_ASSEMBLY_RATIO.out.versions)

        GENERATE_PIPELINE_STATS_WF (
        
            //start changes here
            FASTP_TRIMD.out.reads, \
            GATHERING_READ_QC_STATS.out.fastp_raw_qc, \
            GATHERING_READ_QC_STATS.out.fastp_total_qc, \
            [], \
            KRAKEN2_TRIMD.out.report, \
            KRAKEN2_TRIMD.out.krona_html, \
            KRAKEN2_TRIMD.out.k2_bh_summary, \
            //end changes
            [], \
            [], \
            [], \
            [], \
            [], \
            [], \
            [], \
            RENAME_FASTA_HEADERS.out.renamed_scaffolds, \
            BBMAP_REFORMAT.out.filtered_scaffolds, \
            MLST.out.tsv, \
            GAMMA_HV.out.gamma, \
            GAMMA_AR.out.gamma, \
            GAMMA_PF.out.gamma, \
            QUAST.out.report_tsv, \
            [], [], [], [], \
            KRAKEN2_WTASMBLD.out.report, \
            KRAKEN2_WTASMBLD.out.krona_html, \
            KRAKEN2_WTASMBLD.out.k2_bh_summary, \
            DETERMINE_TAXA_ID.out.taxonomy, \
            FORMAT_ANI.out.ani_best_hit, \
            CALCULATE_ASSEMBLY_RATIO.out.ratio, \
            AMRFINDERPLUS_RUN.out.report, \
            CALCULATE_ASSEMBLY_RATIO.out.gc_content, \
            false
        )

        // Combining output based on meta.id to create summary by sample -- is this verbose, ugly and annoying? yes, if anyone has a slicker way to do this we welcome the input. 
        line_summary_ch = MLST.out.tsv.map{                              meta, tsv             -> [[id:meta.id], tsv]},             by: [0])\
        .join(GAMMA_HV.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},           by: [0])\
        .join(GAMMA_AR.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},           by: [0])\
        .join(GAMMA_PF.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},           by: [0])\
        .join(QUAST.out.report_tsv.map{                                  meta, report_tsv      -> [[id:meta.id], report_tsv]},      by: [0])\
        .join(CALCULATE_ASSEMBLY_RATIO.out.ratio.map{                    meta, ratio           -> [[id:meta.id], ratio]},           by: [0])\
        .join(GENERATE_PIPELINE_STATS_WF.out.pipeline_stats.map{         meta, pipeline_stats  -> [[id:meta.id], pipeline_stats]},  by: [0])\
        .join(DETERMINE_TAXA_ID.out.taxonomy.map{                        meta, taxonomy        -> [[id:meta.id], taxonomy]},        by: [0])\
        //.join(KRAKEN2_TRIMD.out.k2_bh_summary.map{                       meta, k2_bh_summary   -> [[id:meta.id], k2_bh_summary]},   by: [0])\
        .join(AMRFINDERPLUS_RUN.out.report.map{                          meta, report          -> [[id:meta.id], report]}, by: [0])

        // Generate summary per sample that passed SPAdes
        CREATE_SUMMARY_LINE(
            line_summary_ch
        )
        ch_versions = ch_versions.mix(CREATE_SUMMARY_LINE.out.versions)

        // Collect all the summary files prior to fetch step to force the fetch process to wait
        summaries_ch                = CREATE_SUMMARY_LINE.out.line_summary.collect()

        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            all_summaries_ch, false
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

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)

    emit:
        scaffolds        = BBMAP_REFORMAT.out.filtered_scaffolds
        mlst             = MLST.out.tsv
        amrfinder_report = AMRFINDERPLUS_RUN.out.report
        gamma_ar         = GAMMA_AR.out.gamma
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (count == 0){
        if (params.email || params.email_on_fail) {
            NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
        }
        NfcoreTemplate.summary(workflow, params, log)
        count++
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
