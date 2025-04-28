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

include { GET_RAW_STATS         } from '../modules/local/get_raw_stats'
include { CORRUPTION_CHECK      } from '../modules/local/fairy_corruption_check'
include { BBDUK                 } from '../modules/local/bbduk'
include { FASTP as FASTP_LR     } from '../modules/local/fastp'
include { NANOQ                 } from '../modules/local/long_read/nanoq'
include { RASUSA                } from '../modules/local/long_read/rasusa'
include { UNICYCLER             } from '../modules/local/long_read/unicycler'
include { CIRCLATOR             } from '../modules/local/long_read/circlator'
include { BWA                   } from '../modules/local/long_read/bwa'
include { POLYPOLISH            } from '../modules/local/long_read/polypolish'
include { BANDAGE               } from '../modules/local/long_read/bandage'
include { CREATE_SAMPLESHEET    } from '../modules/local/create_samplesheet'
include { FASTQC                } from '../modules/local/fastqc'
include { RAWSTATS              } from '../modules/local/long_read/seqkit'
include { LRGE                  } from '../modules/local/long_read/estimation'


/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/

include { INPUT_CHECK                    } from '../subworkflows/local/input_check'

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

workflow PHOENIX_HYBRID_WF {

    take:
        ch_input

    main:
        ch_versions = Channel.empty() // Used to collect the software versions
        // Allow outdir to be relative
        outdir_path = Channel.fromPath(params.outdir, relative: true)

        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //
        INPUT_CHECK (
            ch_input, "hybrid"
        )
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        //fairy compressed file corruption check & generate read stats
        CORRUPTION_CHECK (
            INPUT_CHECK.out.reads, false // true says busco is being run in this workflow
        )
        ch_versions = ch_versions.mix(CORRUPTION_CHECK.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Short and Long Read trimming with Long Read subsampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

        //Combining reads with output of corruption check. By=2 is for getting R1 and R2 results
        //The mapping here is just to get things in the right bracket so we can call var[0]
        read_stats_ch = INPUT_CHECK.out.reads.join(CORRUPTION_CHECK.out.outcome_to_edit, by: [0,0]) 
            .join(CORRUPTION_CHECK.out.outcome_to_edit.splitCsv(strip:true, by:2).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0]]]}, by: [0,0])
        //.filter{ it[3].findAll {!it.contains('FAILED')}}

        //Get stats on raw reads if the reads aren't corrupted
        GET_RAW_STATS (
            read_stats_ch, false // false says no busco is being run
        )
        ch_versions = ch_versions.mix(GET_RAW_STATS.out.versions)

        // Combining reads with output of corruption check
        bbduk_ch = INPUT_CHECK.out.reads.join(GET_RAW_STATS.out.outcome_to_edit.splitCsv(strip:true, by:3).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0]]]}, by: [0,0])

        // Remove PhiX reads
        BBDUK (
            bbduk_ch, params.bbdukdb
        )
        ch_versions = ch_versions.mix(BBDUK.out.versions)

        // Trim and remove low quality reads
        FASTP_LR (
            BBDUK.out.reads, true, false
        )
        ch_versions = ch_versions.mix(FASTP_LR.out.versions)

        // Long read subsampling and QC
        RAWSTATS(INPUT_CHECK.out.long_read)
        //ch_versions = ch_versions.mix(RAWSTATS.out.versions.first())
    
        LRGE(RAWSTATS.out.fastq_lr)
        //ch_versions = ch_versions.mix(LRGE.out.versions)

        sub_ch = RAWSTATS.out.fastq_lr.map{    meta, fastq_lr       -> [meta, fastq_lr]}\
        .join(LRGE.out.estimation.map{                   meta, estimation            -> [meta, estimation]}, by: [0])
        
        RASUSA (sub_ch,params.depth)
        //ch_versions = ch_versions.mix(RASUSA.out.versions)

        stat_ch = RAWSTATS.out.rawstats.map{    meta, rawstats       -> [meta, rawstats]}\
        .join(RASUSA.out.subfastq.map{                   meta, subfastq            -> [meta, subfastq]}, by: [0])
        
        NANOQ (stat_ch,params.length,params.qscore)
        //ch_versions = ch_versions.mix(NANOQ.out.versions)

    
        /*FASTQC (
            INPUT_CHECK.out.reads
        )
        ch_versions = ch_versions.mix(FASTQC.out.versions)*/

/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    De novo Assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
     
    UNICYCLER (
        NANOQ.out.fastq
    )
    ch_versions = ch_versions.mix(UNICYCLER.out.versions)

/*

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Genome polishing using short reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    CIRCLATOR (
        UNICYCLER.out.fasta
    )
    ch_versions = ch_versions.mix(CIRCLATOR.out.versions)
    
    bwa_ch = CIRCLATOR.out.fasta.map{meta, fasta -> [ meta, fasta ]}\
        .join(FASTP_LR.out.reads.map{meta, reads -> [ [id:meta.id, single_end:true], reads ]}, by: [[0][0],[0][1]])

    BWA (
        bwa_ch
    )
    ch_versions = ch_versions.mix(BWA.out.versions)

    polish_ch = CIRCLATOR.out.fasta.map{meta, fasta -> [ meta, fasta ]}\
        .join(BWA.out.sam.map{          meta, sam   -> [ meta, sam ]}, by: [[0][0],[0][1]])

    POLYPOLISH (
        polish_ch
    )
    ch_versions = ch_versions.mix(POLYPOLISH.out.versions)
   
    BANDAGE (
        UNICYCLER.out.gfa
    )
    ch_versions = ch_versions.mix(BANDAGE.out.versions)

    emit:
        // emits should either be a scaffolds or samplesheet, see comments in main nf.
        scaffolds         = POLYPOLISH.out.assembly.collect()
        valid_samplesheet = INPUT_CHECK.out.valid_samplesheet
        nanostat          = NANOQ.out.nano_stats
        versions          = ch_versions
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

/* Adding if/else for running on ICA
if (params.ica==false) {
    workflow.onComplete {
        if (count == 0) {
            if (params.email || params.email_on_fail) {
                NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
            }
            NfcoreTemplate.summary(workflow, params, log)
            count++
        }
    }
} else if (params.ica==true) {
    workflow.onComplete { 
        if (workflow.success) {
            println("Pipeline Completed Successfully")
            NfcoreTemplate.summary(workflow, params, log)
            System.exit(0)
        } else {
            println("Pipeline Completed with Errors")
            NfcoreTemplate.summary(workflow, params, log)
            System.exit(1)
        }
    }
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