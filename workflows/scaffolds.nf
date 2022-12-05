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

include { ASSET_CHECK                    } from '../modules/local/asset_check'
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
        
        spades_ch = SCAFFOLDS_INPUT_CHECK.out.scaffolds.map{meta, scaffolds -> [ [id:meta.id, single_end:true], scaffolds]}

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
