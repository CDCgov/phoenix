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
include { BBMAP_REFORMAT                 } from '../modules/local/contig_less500'
include { QUAST                          } from '../modules/local/quast'
include { MASH_DIST                      } from '../modules/local/mash_distance'
include { FASTANI                        } from '../modules/local/fastani'
include { DETERMINE_TOP_TAXA             } from '../modules/local/determine_top_taxa'
include { FORMAT_ANI                     } from '../modules/local/format_ANI_best_hit'
include { DETERMINE_TAXA_ID              } from '../modules/local/tax_classifier'
include { PROKKA                         } from '../modules/local/prokka'
include { GET_TAXA_FOR_AMRFINDER         } from '../modules/local/get_taxa_for_amrfinder'
include { AMRFINDERPLUS_RUN              } from '../modules/local/run_amrfinder'
include { CALCULATE_ASSEMBLY_RATIO       } from '../modules/local/assembly_ratio'
include { CREATE_SUMMARY_LINE            } from '../modules/local/phoenix_summary_line'
include { FETCH_FAILED_SUMMARIES         } from '../modules/local/fetch_failed_summaries'
include { GATHER_SUMMARY_LINES           } from '../modules/local/phoenix_summary'
include { GRIPHIN                        } from '../modules/local/griphin'

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/
include { CREATE_INPUT_CHANNEL                     } from '../subworkflows/local/create_input_channel'
include { ADD_EMPTY_CHANNEL                        } from '../subworkflows/local/add_empty_channel'
include { GENERATE_PIPELINE_STATS_WF               } from '../subworkflows/local/generate_pipeline_stats'
include { KRAKEN2_WF as KRAKEN2_ASMBLD             } from '../subworkflows/local/kraken2krona'
include { KRAKEN2_WF as KRAKEN2_WTASMBLD           } from '../subworkflows/local/kraken2krona'
include { DO_MLST                                  } from '../subworkflows/local/do_mlst'

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
    GROOVY FUNCTIONS
========================================================================================
*/

def add_empty_ch(input_ch) {
    meta_id = input_ch[0]
    output_array = [ meta_id, input_ch[1], input_ch[2], []]
    return output_array
}

// Groovy funtion to make [ meta.id, [] ] - just an empty channel
def create_empty_ch(input_for_meta) { // We need meta.id associated with the empty list which is why .ifempty([]) won't work
    meta_id = input_for_meta[0]
    output_array = [ meta_id, [] ]
    return output_array
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow SCAFFOLDS_EXQC {
    take:
        ch_input
        ch_input_indir

    main:
        ch_versions = Channel.empty() // Used to collect the software versions
        // Allow outdir to be relative
        outdir_path = Channel.fromPath(params.outdir, relative: true)

        CREATE_INPUT_CHANNEL (
            ch_input_indir, ch_input
        )

        //unzip any zipped databases
        ASSET_CHECK (
            params.zipped_sketch
        )

        // Rename scaffold headers
        RENAME_FASTA_HEADERS (
            CREATE_INPUT_CHANNEL.out.scaffolds_ch
        )
        ch_versions = ch_versions.mix(RENAME_FASTA_HEADERS.out.versions)

        // Removing scaffolds <500bp
        BBMAP_REFORMAT (
            RENAME_FASTA_HEADERS.out.renamed_scaffolds
        )
        ch_versions = ch_versions.mix(BBMAP_REFORMAT.out.versions)

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

        if (params.busco_db_path != null) {
            // Allow relative paths for krakendb argument
            busco_db_path = Channel.fromPath(params.busco_db_path, relative: true) 
            busco_ch = BBMAP_REFORMAT.out.filtered_scaffolds.combine(busco_db_path) // Add in krakendb into the fasta channel so each fasta has a krakendb to go with it. 
        } else {
            // passing empty channel for busco db to align with expected inputs for the module
            busco_ch = BBMAP_REFORMAT.out.filtered_scaffolds.map{ meta, scaffolds -> [ [id:meta.id, single_end:meta.single_end], scaffolds, []]}
        }

        // Checking single copy genes for assembly completeness
        BUSCO (
            busco_ch, 'auto', []
        )
        ch_versions = ch_versions.mix(BUSCO.out.versions)

        // Checking for Contamination in assembly creating krona plots and best hit files
        KRAKEN2_ASMBLD (
            BBMAP_REFORMAT.out.filtered_scaffolds,"asmbld", [], QUAST.out.report_tsv
        )
        ch_versions = ch_versions.mix(KRAKEN2_ASMBLD.out.versions)

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
        top_taxa_ch = MASH_DIST.out.dist.map{           meta, dist  -> [[id:meta.id], dist]}\
        .join(BBMAP_REFORMAT.out.filtered_scaffolds.map{meta, reads -> [[id:meta.id], reads ]}, by: [0])

        // Generate file with list of paths of top taxa for fastANI
        DETERMINE_TOP_TAXA (
            top_taxa_ch
        )
        ch_versions = ch_versions.mix(DETERMINE_TOP_TAXA.out.versions)

        // Combining filtered scaffolds with the top taxa list based on meta.id
        top_taxa_list_ch = BBMAP_REFORMAT.out.filtered_scaffolds.map{meta, reads           -> [[id:meta.id], reads]}\
        .join(DETERMINE_TOP_TAXA.out.top_taxa_list.map{              meta, top_taxa_list   -> [[id:meta.id], top_taxa_list ]},   by: [0])\
        .join(DETERMINE_TOP_TAXA.out.reference_files.map{            meta, reference_files -> [[id:meta.id], reference_files ]}, by: [0])

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
        .join(FORMAT_ANI.out.ani_best_hit.map{        meta, ani_best_hit           -> [[id:meta.id], ani_best_hit ]},  by: [0]).map{ it -> add_empty_ch(it) }

        // Getting ID from either FastANI or if fails, from Kraken2
        DETERMINE_TAXA_ID (
            best_hit_ch, params.taxa
        )
        ch_versions = ch_versions.mix(DETERMINE_TAXA_ID.out.versions)
 
        // Perform MLST steps on isolates (with srst2 on internal samples)
        DO_MLST (
            BBMAP_REFORMAT.out.filtered_scaffolds, \
            [], \
            DETERMINE_TAXA_ID.out.taxonomy, \
            false
        )

        // get gff and protein files for amrfinder+
        PROKKA (
            BBMAP_REFORMAT.out.filtered_scaffolds, [], []
        )
        ch_versions = ch_versions.mix(PROKKA.out.versions)

        /*/ Fetch AMRFinder Database
        AMRFINDERPLUS_UPDATE( )
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)*/

        // Create file that has the organism name to pass to AMRFinder
        GET_TAXA_FOR_AMRFINDER (
            DETERMINE_TAXA_ID.out.taxonomy
        )

        // Combining taxa and scaffolds to run amrfinder and get the point mutations.
        amr_channel = BBMAP_REFORMAT.out.filtered_scaffolds.map{                 meta, reads          -> [[id:meta.id], reads]}\
        .join(GET_TAXA_FOR_AMRFINDER.out.amrfinder_taxa.splitCsv(strip:true).map{meta, amrfinder_taxa -> [[id:meta.id], amrfinder_taxa ]}, by: [0])\
        .join(PROKKA.out.faa.map{                                                meta, faa            -> [[id:meta.id], faa ]},            by: [0])\
        .join(PROKKA.out.gff.map{                                                meta, gff            -> [[id:meta.id], gff ]},            by: [0])

        // Run AMRFinder
        AMRFINDERPLUS_RUN (
            amr_channel, params.amrfinder_db
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
            [], \
            [], \
            [], \
            [], \
            [], \
            [], \
            RENAME_FASTA_HEADERS.out.renamed_scaffolds, \
            BBMAP_REFORMAT.out.filtered_scaffolds, \
            DO_MLST.out.checked_MLSTs, \
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
            AMRFINDERPLUS_RUN.out.mutation_report, \
            CALCULATE_ASSEMBLY_RATIO.out.gc_content, \
            false
        )

        // Creating empty channel that has the form [ meta.id, [] ] that can be passed as a blank below
        empty_ch = RENAME_FASTA_HEADERS.out.renamed_scaffolds.map{ it -> create_empty_ch(it) }

        // Combining output based on meta.id to create summary by sample -- is this verbose, ugly and annoying? yes, if anyone has a slicker way to do this we welcome the input.
        line_summary_ch = empty_ch.map{                                  meta, list            -> [[id:meta.id], list]}\
        .join(DO_MLST.out.checked_MLSTs.map{                             meta, checked_MLSTs   -> [[id:meta.id], checked_MLSTs]},  by: [0])\
        .join(GAMMA_HV.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},          by: [0])\
        .join(GAMMA_AR.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},          by: [0])\
        .join(GAMMA_PF.out.gamma.map{                                    meta, gamma           -> [[id:meta.id], gamma]},          by: [0])\
        .join(QUAST.out.report_tsv.map{                                  meta, report_tsv      -> [[id:meta.id], report_tsv]},     by: [0])\
        .join(CALCULATE_ASSEMBLY_RATIO.out.ratio.map{                    meta, ratio           -> [[id:meta.id], ratio]},          by: [0])\
        .join(GENERATE_PIPELINE_STATS_WF.out.pipeline_stats.map{         meta, pipeline_stats  -> [[id:meta.id], pipeline_stats]}, by: [0])\
        .join(DETERMINE_TAXA_ID.out.taxonomy.map{                        meta, taxonomy        -> [[id:meta.id], taxonomy]},       by: [0])\
        .join(empty_ch.map{                                              meta, list            -> [[id:meta.id], list]},           by: [0])\
        .join(AMRFINDERPLUS_RUN.out.report.map{                          meta, report          -> [[id:meta.id], report]},         by: [0])

        // Generate summary per sample that passed SPAdes
        CREATE_SUMMARY_LINE(
            line_summary_ch
        )
        ch_versions = ch_versions.mix(CREATE_SUMMARY_LINE.out.versions)

        // Collect all the summary files prior to fetch step to force the fetch process to wait
         summaries_ch = CREATE_SUMMARY_LINE.out.line_summary.collect()

        // Combining sample summaries into final report
        GATHER_SUMMARY_LINES (
            summaries_ch, outdir_path, false
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        GRIPHIN (
            summaries_ch, CREATE_INPUT_CHANNEL.out.valid_samplesheet, params.ardb, outdir_path
        )
        ch_versions = ch_versions.mix(GATHER_SUMMARY_LINES.out.versions)

        // Collecting the software versions
        CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

    emit:
        scaffolds        = BBMAP_REFORMAT.out.filtered_scaffolds
        mlst             = DO_MLST.out.checked_MLSTs
        amrfinder_report = AMRFINDERPLUS_RUN.out.report
        gamma_ar         = GAMMA_AR.out.gamma
        summary_report   = GATHER_SUMMARY_LINES.out.summary_report
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
