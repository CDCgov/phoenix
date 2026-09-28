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
//def multiqc_report = []

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

include { SEQKIT_RAWSTATS       } from '../modules/local/long_read/seqkit'
include { CORRUPTION_CHECK      } from '../modules/local/fairy_corruption_check'
include { NANOQ                 } from '../modules/local/long_read/nanoq'
include { RASUSA                } from '../modules/local/long_read/rasusa'
include { FLYE                  } from '../modules/local/long_read/flye'
include { MEDAKA                } from '../modules/local/long_read/medaka'
include { CIRCLATOR             } from '../modules/local/long_read/circlator'
include { BANDAGE               } from '../modules/local/long_read/bandage'
include { CREATE_SAMPLESHEET    } from '../modules/local/create_samplesheet'
include { LRGE                  } from '../modules/local/long_read/estimation'
include { PLSDB_ASSET_CHECK     } from '../modules/local/long_read/plsdb_asset_check'
include { FASTQC                } from '../modules/local/fastqc'

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK              } from '../subworkflows/local/input_check'
include { PLASMID_CHARACTERIZATION } from '../subworkflows/local/plasmid_characterization'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//

//include { MULTIQC                      } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow PHOENIX_LR_WF {

    take:
        ch_input

    main:
        ch_versions = Channel.empty() // Used to collect the software versions
        // Allow outdir to be relative
        outdir_path = Channel.fromPath(params.outdir, relative: true)

        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        INPUT_CHECK (
            ch_input, false
        )
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        // Handling plsdb database for plasmid characterization subworkflow - check if plsdb exists; if not, download and makeblastdb
        PLSDB_ASSET_CHECK (
            params.plsdb_dir
        )
        ch_versions = ch_versions.mix(PLSDB_ASSET_CHECK.out.versions)

        //fairy compressed file corruption check & generate read stats
        CORRUPTION_CHECK (
            INPUT_CHECK.out.reads, false, workflow.manifest.version  // true says busco is being run in this workflow
        )
        ch_versions = ch_versions.mix(CORRUPTION_CHECK.out.versions)

        raw_reads_ch = INPUT_CHECK.out.reads.join(CORRUPTION_CHECK.out.outcome_to_edit, by: [0,0])
                        .join(CORRUPTION_CHECK.out.outcome_to_edit.splitCsv(strip:true, by:2).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0]]]}, by: [0,0])
                        .filter { meta, reads, fairy_outcome_to_edit, fairy_outcome -> fairy_outcome.every { it.startsWith("PASSED:") } } //if the files are not corrupt then get the read stats
                        .map{ meta, reads, fairy_outcome_to_edit, fairy_outcome -> return [meta, reads, fairy_outcome_to_edit] }

        SEQKIT_RAWSTATS (
            raw_reads_ch
        )
        ch_versions = ch_versions.mix(SEQKIT_RAWSTATS.out.versions)

        raw_reads_ch2 = INPUT_CHECK.out.reads.join(SEQKIT_RAWSTATS.out.outcome_to_edit, by: [0,0])
                        .join(SEQKIT_RAWSTATS.out.outcome_to_edit.splitCsv(strip:true, by:2).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0]]]}, by: [0,0])
                        .filter { meta, reads, fairy_outcome_to_edit, fairy_outcome -> fairy_outcome.every { it.startsWith("PASSED:") } } //if the files are not corrupt then get the read stats
                        .map{ meta, reads, fairy_outcome_to_edit, fairy_outcome -> return [meta, reads, fairy_outcome_to_edit] }

        // Estimating genome size from long read overlaps
        LRGE (
            raw_reads_ch2
        )
        ch_versions = ch_versions.mix(LRGE.out.versions)

        sub_ch = INPUT_CHECK.out.reads.map{meta, reads -> [meta, reads]}.join(LRGE.out.estimation.map{meta, estimation -> [meta, estimation]}, by: [0,0]).join(LRGE.out.outcome_to_edit, by: [0,0])
                    .join(LRGE.out.outcome_to_edit.splitCsv(strip:true, by:3).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0]]]}, by: [0,0])
                    .filter{ meta, reads, estimation, fairy_outcome_to_edit, fairy_outcome -> fairy_outcome.every { it.startsWith("PASSED:") } } //if the files are not corrupt then get the read stats
                    .map{ meta, reads, estimation, fairy_outcome_to_edit, fairy_outcome -> return [meta, reads, estimation] }

        //subsample reads to a target depth of coverage using RASUSA
        RASUSA (
            sub_ch,
            params.depth
        )
        ch_versions = ch_versions.mix(RASUSA.out.versions)

        // RAWSTATS.out.rawstats -
        //The awk command reads that file, goes to line 2 (NR==2), and prints the 4th whitespace-separated field ($4). This is a value extracted from raw-read statistics
        /*echo -e "raw_reads trim_reads bases n50 longest shortest mean_length median_length mean_quality median_quality\n\$(awk 'NR==2 {print \$4}' $rawstats) \$(awk 'NR==2 {print}' trim.txt)" > ${meta.id}_stats_mqc.txt
        # convert to stats for reporting in griphin later to csv
        sed 's/ /,/g' ${meta.id}_stats_mqc.txt > ${meta.id}_nanoq_stats.csv */

        stat_ch = RASUSA.out.subfastq.map{meta, subfastq -> [meta, subfastq]}.join(LRGE.out.outcome_to_edit, by: [0,0])
                    .join(LRGE.out.outcome_to_edit.splitCsv(strip:true, by:3).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0]]]}, by: [0,0])
                    .filter{ meta, subfastq, fairy_outcome_to_edit, fairy_outcome -> fairy_outcome.every { it.startsWith("PASSED:") } } //if the files are not corrupt then get the read stats
                    .map{ meta, subfastq, fairy_outcome_to_edit, fairy_outcome -> return [meta, subfastq, fairy_outcome_to_edit] }

        //Ultra-fast quality control and summary reports for nanopore reads
        NANOQ (
            stat_ch,
            params.length,
            params.qscore
        )
        ch_versions = ch_versions.mix(NANOQ.out.versions)

        trimed_ch = NANOQ.out.fastq.map{meta, fastq -> [meta, fastq]}.join(NANOQ.out.outcome_to_edit, by: [0,0])
                    .join(NANOQ.out.outcome_to_edit.splitCsv(strip:true, by:4).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0]]]}, by: [0,0])
                    .filter{ meta, subfastq, fairy_outcome_to_edit, fairy_outcome -> fairy_outcome.every { it.startsWith("PASSED:") } } //if the files are not corrupt then get the read stats
                    .map{ meta, subfastq, fairy_outcome_to_edit, fairy_outcome -> return [meta, subfastq, fairy_outcome_to_edit] }

        // De novo assembly with Flye
        FLYE (
            trimed_ch
        )
        ch_versions = ch_versions.mix(FLYE.out.versions)

        FLYE.out.outcome.view { meta, file -> "outcome for fly check ${meta.id}: ${file.text}" }

        assembly_to_polished_ch = FLYE.out.fasta.map{meta, fasta -> [meta, fasta]}.join(NANOQ.out.fastq.map{meta, fastq -> [meta, fastq]}, by: [0,0])

        // Medaka long-read polishing
        MEDAKA (
            assembly_to_polished_ch
        )
        ch_versions = ch_versions.mix(MEDAKA.out.versions)

        // Circlator for circularization and orientation
        CIRCLATOR (
            MEDAKA.out.fasta
        )
        ch_versions = ch_versions.mix(CIRCLATOR.out.versions)

        // create bandage graph for visualization of assembly
        BANDAGE (
            FLYE.out.gfa
        )
        ch_versions = ch_versions.mix(BANDAGE.out.versions)

        // Plasmid characterization subworkflow
        PLASMID_CHARACTERIZATION (
            MEDAKA.out.fasta, \
            //SCAFFOLD_COUNT_CHECK.out.outcome, \
            params.plsdb_dir, \
            params.conf, \
            params.viz, \
            //params.plsdbfasta, \
            PLSDB_ASSET_CHECK.out.plsdb, \
            params.vfdb
        )
        ch_versions = ch_versions.mix(PLASMID_CHARACTERIZATION.out.versions)    

    emit:
        // emits should either be a scaffolds or samplesheet, see comments in main nf.
        scaffolds             = MEDAKA.out.fasta_gz.collect()
        valid_samplesheet     = INPUT_CHECK.out.valid_samplesheet
        trimmed_stats         = NANOQ.out.trimmed_stats
        raw_stats             = SEQKIT_RAWSTATS.out.rawstats
        fairy_outcome_to_edit = FLYE.out.outcome
        versions              = ch_versions
        // emits for plasmid characterization outputs
        // plasmidID        = PLASMID_CHARACTERIZATION.out.blast_out
        // plasmidANI       = PLASMID_CHARACTERIZATION.out.ani_out
        // plasmidVF        = PLASMID_CHARACTERIZATION.out.plasmid_vf_out

}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

// Adding if/else for running on ICA
if (params.ica==false) {
    // do nothing, not running ICA and no erros occurred
} else if (params.ica==true) {
    workflow.onError { 
        // copy intermediate files + directories
        println("Getting intermediate files from ICA")
        ['cp','-r',"${workflow.workDir}","${workflow.launchDir}/out"].execute()
        // return trace files
        println("Returning workflow run-metric reports from ICA")
        ['find','/ces','-type','f','-name','\"*.ica\"','2>','/dev/null', '|', 'grep','"report"' ,'|','xargs','-i','cp','-r','{}',"${workflow.launchDir}/out"].execute()
    }
} else {
        error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods."
}

/*
========================================================================================
    THE END
========================================================================================

}
*/
