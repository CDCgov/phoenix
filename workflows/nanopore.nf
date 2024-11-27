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

include { SAMPLESHEET_CHECK     } from '../modules/local/long_read/check_lr' // make this folder and add module files fill xxx with the name of the file
include { NANOQ                 } from '../modules/local/long_read/nanoq'
include { RASUSA                } from '../modules/local/long_read/rasusa'
include { FLYE                  } from '../modules/local/long_read/flye'
include { MEDAKA                } from '../modules/local/long_read/medaka'
include { CIRCLATOR             } from '../modules/local/long_read/circlator'
include { BANDAGE               } from '../modules/local/long_read/bandage'
//include { PLASMID                 } from '../modules/local/long_read/plasmid'
//include { PLASMID_AR            } from '../modules/local/long_read/plasmid_ar'
//include { PLASMID_VF            } from '../modules/local/long_read/plasmid_vf'


//include { FLYE            } from '../modules/local/long_read/xxx'
//include { CIRCLATOR       } from '../modules/local/long_read/xxx'
// fill out the rest ... 

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/

include { LR_CHECK           } from '../subworkflows/local/lr_check'

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
        // Allow relative paths for krakendb argument
        kraken2_db_path  = Channel.fromPath(params.kraken2db, relative: true)

        //
        // SUBWORKFLOW: Read in samplesheet, validate and stage input files
        //
        LR_CHECK (
            file(ch_input)
        )
        ch_versions = ch_versions.mix(LR_CHECK.out.versions)
    
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Read Preprocessing
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */  
    NANOQ (LR_CHECK.out.reads)
    ch_versions = ch_versions.mix(NANOQ.out.versions.first())

    //FASTQC (LR_CHECK.out.reads)
    //ch_versions = ch_versions.mix(FASTQC.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Read Subsampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    RASUSA (NANOQ.out.fastq)
    ch_versions = ch_versions.mix(RASUSA.out.versions)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    De novo Assembly
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
    FLYE (RASUSA.out.fastq)
    ch_versions = ch_versions.mix(FLYE.out.versions)
   
    CIRCLATOR (FLYE.out.fasta)
    ch_versions = ch_versions.mix(CIRCLATOR.out.versions)

    BANDAGE (FLYE.out.gfa)
    ch_versions = ch_versions.mix(BANDAGE.out.versions)

    polish_ch = CIRCLATOR.out.fasta.map       {meta,fasta                 -> [meta,fasta]}\
        .join(NANOQ.out.fastq.map{   meta, fastq  -> [meta, fastq]},    by: [0])\


    MEDAKA (polish_ch)
    ch_versions = ch_versions.mix(MEDAKA.out.versions)
    

    /*plasmid_ch = PLASMIDmarker .out.tsv.map  {meta,tsv                -> [meta,tsv]}\
        .join(MEDAKA.out.fasta.map{   meta, fasta  -> [meta, fasta]},    by: [0])\


    //This is preliminary and in development
    markerANI (plasmid_ch,params.viz)*/
/*
    PLASMID (MEDAKA.out.fasta,params.plsdb,params.conf) 
    PLASMID_AR (PLASMID.out.contigs,params.ardb)
    PLASMID_VF (PLASMID.out.contigs,params.hvdb)
*/
    emit:
        // emits should either be a scaffolds or samplesheet, see comments in main nf.
        scaffolds = MEDAKA.out.fasta_fin // or whatever the best scaffold file is... 


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