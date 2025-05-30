#!/usr/bin/env nextflow
/*
========================================================================================
    CDCgov/phoenix
========================================================================================
    Github : https://github.com/CDCgov/phoenix
    Slack  : https://staph-b-dev.slack.com/channels/phoenix-dev
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// ANSI escape code for orange (bright yellow)
def orange = '\033[38;5;208m'
def reset = '\033[0m'

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

//Check coverage is above its threshold
if (params.coverage.toInteger() < 30) { exit 1, 'The minimum coverage allowed for QA/QC purposes is 30 and is the default. Please choose a value >=30.' }
// Check for incorrect --output parameter
params.output = "" /// Initialise param so no warning is printed
if (params.output) { exit 1, "ERROR: Unknown parameter '--output'. Did you mean '--outdir'?" }
//comment out in v2.3.0 to run --centar
if (params.centar == true) { exit 1, "Sorry, --centar available yet as it's validation isn't complete. It will be released with a newer version of phx in the future." }

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { PHOENIX_EXTERNAL            } from './workflows/phoenix'
include { PHOENIX_EXQC                } from './workflows/cdc_phoenix'
include { SCAFFOLDS_EXTERNAL          } from './workflows/scaffolds'
include { SCAFFOLDS_EXQC              } from './workflows/cdc_scaffolds'
include { SRA_PREP                    } from './workflows/sra_prep'
include { CLIA_INTERNAL               } from './workflows/clia'
include { UPDATE_PHOENIX_WF           } from './workflows/update_phoenix'
include { RUN_CENTAR                  } from './workflows/centar'
include { COMBINE_GRIPHINS_WF         } from './workflows/combine_griphins'

//
// WORKFLOW: Run main cdcgov/phoenix analysis pipeline
//
workflow PHOENIX {
    //Check path of kraken2db
    if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified! Use --kraken2db to tell PHoeNIx where to find the database.' }

    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db] //removed , params.fasta to stop issue w/connecting to aws and igenomes not used
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters

    //input on command line
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'For -entry PHOENIX: Input samplesheet not specified!' }
    ch_versions = Channel.empty() // Used to collect the software versions

    main:
        PHOENIX_EXTERNAL ( ch_input, ch_versions, true, params.centar )

    emit:
        scaffolds        = PHOENIX_EXTERNAL.out.scaffolds
        trimmed_reads    = PHOENIX_EXTERNAL.out.trimmed_reads
        mlst             = PHOENIX_EXTERNAL.out.mlst
        amrfinder_output = PHOENIX_EXTERNAL.out.amrfinder_output
        gamma_ar         = PHOENIX_EXTERNAL.out.gamma_ar
        phx_summary      = PHOENIX_EXTERNAL.out.phx_summary
        //output for phylophoenix
        griphin_tsv      = PHOENIX_EXTERNAL.out.griphin_tsv
        griphin_excel    = PHOENIX_EXTERNAL.out.griphin_excel
        dir_samplesheet  = PHOENIX_EXTERNAL.out.dir_samplesheet
        //output for ncbi upload 
        ncbi_sra_sheet       = params.create_ncbi_sheet ? PHOENIX_EXTERNAL.out.ncbi_sra_sheet : null
        ncbi_biosample_sheet = params.create_ncbi_sheet ? PHOENIX_EXTERNAL.out.ncbi_biosample_sheet : null
}

//
// WORKFLOW: Run internal version of cdcgov/phoenix analysis pipeline that includes BUSCO, SRST2 and KRAKEN_ASMBLED
//
workflow CDC_PHOENIX {
    //Check path of kraken2db
    if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified! Use --kraken2db to tell PHoeNIx where to find the database.' }

    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters

    //input on command line
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'For -entry CDC_PHOENIX: Input samplesheet not specified!' }
    ch_versions = Channel.empty() // Used to collect the software versions

    main:
        PHOENIX_EXQC ( ch_input, ch_versions, true, params.centar )

    emit:
        scaffolds        = PHOENIX_EXQC.out.scaffolds
        trimmed_reads    = PHOENIX_EXQC.out.trimmed_reads
        mlst             = PHOENIX_EXQC.out.mlst
        amrfinder_output = PHOENIX_EXQC.out.amrfinder_output
        gamma_ar         = PHOENIX_EXQC.out.gamma_ar
        phx_summary      = PHOENIX_EXQC.out.phx_summary
        //output for phylophoenix
        griphin_tsv      = PHOENIX_EXQC.out.griphin_tsv
        griphin_excel    = PHOENIX_EXQC.out.griphin_excel
        dir_samplesheet  = PHOENIX_EXQC.out.dir_samplesheet
        //output for ncbi upload 
        ncbi_sra_sheet       = params.create_ncbi_sheet ? PHOENIX_EXQC.out.ncbi_sra_sheet : null
        ncbi_biosample_sheet = params.create_ncbi_sheet ? PHOENIX_EXQC.out.ncbi_biosample_sheet : null
}

/*
========================================================================================
    RUN SRA WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Run internal version of phx based on sample SRA names, these will be pulled from NCBI for you. 
//
workflow SRA {
    //Check path of kraken2db
    if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified! Use --kraken2db to tell PHoeNIx where to find the database.' }

    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input_sra, params.multiqc_config, params.kraken2db ]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Checking that --create_ncbi_sheet wasn't passed
    if (params.create_ncbi_sheet) { exit 1, '--create_ncbi_sheet is not a valid argument for -entry SRA.' }

    // Check mandatory parameters
    //input on command line
    if (params.input_sra) {
        //create channel input
        ch_input = file(params.input_sra)
        //Check that SRR numbers are passed not SRX
        if (ch_input) {
            // Read the contents of the file
            def sraNumbers = ch_input.text.readLines()
            // Check each line in the file
            for (sraNumber in sraNumbers) {
                // Check if it starts with "SRR"
                if (!sraNumber.startsWith("SRR")) {
                    exit 1, "Invalid value in ${params.input_sra}. Only SRR numbers are allowed for -entry SRA, but found: $sraNumber"
                }
            }
        }
    } else { exit 1, 'For -entry SRA: Input samplesheet not specified! Make sure to use --input_sra NOT --input' }

    main:
        // pull data and create samplesheet for it.
        SRA_PREP ( ch_input )
        // pass samplesheet to PHOENIX
        PHOENIX_EXTERNAL ( SRA_PREP.out.samplesheet, SRA_PREP.out.versions, false, params.centar )

    emit:
        scaffolds        = PHOENIX_EXTERNAL.out.scaffolds
        trimmed_reads    = PHOENIX_EXTERNAL.out.trimmed_reads
        mlst             = PHOENIX_EXTERNAL.out.mlst
        amrfinder_output = PHOENIX_EXTERNAL.out.amrfinder_output
        gamma_ar         = PHOENIX_EXTERNAL.out.gamma_ar
        phx_summary      = PHOENIX_EXTERNAL.out.phx_summary
        //output for phylophoenix
        griphin_tsv      = PHOENIX_EXTERNAL.out.griphin_tsv
        griphin_excel    = PHOENIX_EXTERNAL.out.griphin_excel
        dir_samplesheet  = PHOENIX_EXTERNAL.out.dir_samplesheet
}

//
// WORKFLOW: Run cdc version of phx based on sample SRA names, the fastq will be pulled from NCBI for you. 
//
workflow CDC_SRA {
    //Check path of kraken2db
    if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified! Use --kraken2db to tell PHoeNIx where to find the database.' }

    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input_sra, params.multiqc_config, params.kraken2db]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Checking that --create_ncbi_sheet wasn't passed
    if (params.create_ncbi_sheet) { exit 1, '--create_ncbi_sheet is not a valid argument for -entry CDC_SRA.' }

    // Check mandatory parameters
    //input on command line
    if (params.input_sra) {
        //create channel input
        ch_input = file(params.input_sra)
        //Check that SRR numbers are passed not SRX
        if (ch_input) {
            // Read the contents of the file
            def sraNumbers = ch_input.text.readLines()
            // Check each line in the file
            for (sraNumber in sraNumbers) {
                // Check if it starts with "SRR"
                if (!sraNumber.startsWith("SRR")) {
                    exit 1, "Invalid value in ${params.input_sra}. Only SRR numbers are allowed for -entry CDC_SRA, but found: $sraNumber"
                }
            }
        }
    } else { exit 1, 'For -entry CDC_SRA: Input samplesheet not specified! Make sure to use --input_sra NOT --input' }

    main:
        // pull data and create samplesheet for it.
        SRA_PREP ( ch_input )
        // pass samplesheet to PHOENIX
        PHOENIX_EXQC ( SRA_PREP.out.samplesheet, SRA_PREP.out.versions, false, params.centar )

    emit:
        scaffolds        = PHOENIX_EXQC.out.scaffolds
        trimmed_reads    = PHOENIX_EXQC.out.trimmed_reads
        mlst             = PHOENIX_EXQC.out.mlst
        amrfinder_output = PHOENIX_EXQC.out.amrfinder_output
        gamma_ar         = PHOENIX_EXQC.out.gamma_ar
        phx_summary      = PHOENIX_EXQC.out.phx_summary
        //output for phylophoenix
        griphin_tsv      = PHOENIX_EXQC.out.griphin_tsv
        griphin_excel    = PHOENIX_EXQC.out.griphin_excel
        dir_samplesheet  = PHOENIX_EXQC.out.dir_samplesheet
}

/*
========================================================================================
    RUN SCAFFOLD WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Entry point to analyze scaffold file(s) and run everything after Spades
//
workflow SCAFFOLDS {
    //Check path of kraken2db
    if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified! Use --kraken2db to tell PHoeNIx where to find the database.' }

    // Checking that --create_ncbi_sheet wasn't passed
    if (params.create_ncbi_sheet) { exit 1, '--create_ncbi_sheet is not a valid argument for -entry SCAFFOLDS.' }

    // Validate input parameters
    // Check input path parameters to see if they exist
    if (params.input != null ) {  // if a samplesheet is passed
        if (params.indir != null ) { //if samplesheet is passed and an input directory exit
            exit 1, 'For -entry SCAFFOLDS: You need EITHER an input samplesheet or a directory! Just pick one.' 
        } else { // if only samplesheet is passed check to make sure input is an actual file
            def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = null //keep input directory null if not passed
            // get full path for input and make channel
            if (params.input) { ch_input = file(params.input) }
        }
    } else {
        if (params.indir != null ) { // if no samplesheet is passed, but an input directory is given
            ch_input = null //keep samplesheet input null if not passed
            def checkPathParamList = [ params.indir, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = Channel.fromPath(params.indir, relative: true)
        } else { // if no samplesheet is passed and no input directory is given
            exit 1, 'For -entry SCAFFOLDS: You need EITHER an input samplesheet or a directory!' 
        }
    }

    main:
        SCAFFOLDS_EXTERNAL ( ch_input, ch_input_indir, params.centar )

    emit:
        scaffolds        = SCAFFOLDS_EXTERNAL.out.scaffolds
        mlst             = SCAFFOLDS_EXTERNAL.out.mlst
        amrfinder_output = SCAFFOLDS_EXTERNAL.out.amrfinder_output
        gamma_ar         = SCAFFOLDS_EXTERNAL.out.gamma_ar
        phx_summary      = SCAFFOLDS_EXTERNAL.out.phx_summary
}

//
// WORKFLOW: Entry point to analyze scaffold file(s) and run everything after Spades
//
workflow CDC_SCAFFOLDS {
    //Check path of kraken2db
    if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified! Use --kraken2db to tell PHoeNIx where to find the database.' }

    // Checking that --create_ncbi_sheet wasn't passed
    if (params.create_ncbi_sheet) { exit 1, '--create_ncbi_sheet is not a valid argument for -entry CDC_SCAFFOLDS.' }

    // Validate input parameters
    // Check input path parameters to see if they exist
    if (params.input != null ) {  // if a samplesheet is passed
        // allow input to be relative
        //input_samplesheet_path = Channel.fromPath(params.input, relative: true)
        if (params.indir != null ) { //if samplesheet is passed and an input directory exit
            exit 1, 'For -entry CDC_SCAFFOLDS: You need EITHER an input samplesheet or a directory! Just pick one.' 
        } else { // if only samplesheet is passed check to make sure input is an actual file
            def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = null //keep input directory null if not passed
            // get full path for input and make channel
            if (params.input) { ch_input = file(params.input) }
        }
    } else {
        if (params.indir != null ) { // if no samplesheet is passed, but an input directory is given
            ch_input = null //keep samplesheet input null if not passed
            def checkPathParamList = [ params.indir, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = Channel.fromPath(params.indir, relative: true)
        } else { // if no samplesheet is passed and no input directory is given
            exit 1, 'For -entry CDC_SCAFFOLDS: You need EITHER an input samplesheet or a directory!' 
        }
    }

    main:
        SCAFFOLDS_EXQC ( ch_input, ch_input_indir, params.centar )

    emit:
        scaffolds        = SCAFFOLDS_EXQC.out.scaffolds
        mlst             = SCAFFOLDS_EXQC.out.mlst
        amrfinder_output = SCAFFOLDS_EXQC.out.amrfinder_output
        gamma_ar         = SCAFFOLDS_EXQC.out.gamma_ar
        phx_summary      = SCAFFOLDS_EXQC.out.phx_summary
}

/*/
// WORKFLOW: Entry point for CLIA analysis
//
workflow CLIA {
    //Check path of kraken2db
    if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified! Use --kraken2db to tell PHoeNIx where to find the database.' }

    // Checking that --create_ncbi_sheet wasn't passed
    if (params.create_ncbi_sheet) { exit 1, '--create_ncbi_sheet is not a valid argument for -entry CLIA.' }

    //Check that SRR numbers are passed no SRX
    if (params.create_ncbi_sheet) {
        // Read the contents of the file
        def sraNumbers = file(params.create_ncbi_sheet).text.readLines()

        // Check each line in the file
        for (sraNumber in sraNumbers) {
            // Check if it starts with "SRR"
            if (!sraNumber.startsWith("SRR")) {
                exit 1, "Invalid value in ${params.create_ncbi_sheet}. Only SRR numbers are allowed, but found: $sraNumber"
            }
        }
    }

    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters
    //input on command line
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'For -entry CLIA: Input samplesheet not specified!' }
    ch_versions = Channel.empty() // Used to collect the software versions

    // Check that a busco_db_path is passed
    // ; means do nothing as that is correct
    if (params.busco_db_path != null) { ; } else { exit 1, 'For -entry CLIA, BUSCO offline mode is not allowed, please pass a path to --busco_db_path!' }

    main:
        CLIA_INTERNAL ( ch_input, ch_versions )

    emit:
        scaffolds        = CLIA_INTERNAL.out.scaffolds
        trimmed_reads    = CLIA_INTERNAL.out.trimmed_reads
        amrfinder_report = CLIA_INTERNAL.out.amrfinder_report
        summary_report   = CLIA_INTERNAL.out.summary_report
}*/

/*
========================================================================================
    RUN PHX Utilities WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Entry point for updating phoenix mlst and ar output
//
workflow UPDATE_PHOENIX {
    //Check path of kraken2db
    if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified! Use --kraken2db to tell PHoeNIx where to find the database.' }

    //Regardless of what is passed outdir needs to be the same as the input dir 
    //if you don't pass outdir then the indir
    //if (params.outdir == "${launchDir}/phx_output" ) { params.outdir = params.indir } else { println("You didn't specify an outdir so phx assumes its the same as the indir.") }

    // check config file
    if (!workflow.configFiles) {
        error "The -c parameter (configuration file) is missing. If you don't pass this then the default databases for this version of phoenix will be used."
    }

    // Check mandatory parameters
    ch_versions = Channel.empty() // Used to collect the software versions
    // Check input path parameters to see if they exist
    if (params.input != null ) {  // if a samplesheet is passed
        //input_samplesheet_path = Channel.fromPath(params.input, relative: true)
        if (params.indir != null ) { //if samplesheet is passed and an input directory exit
            exit 1, 'For -entry UPDATE_CDC_PHOENIX: You need EITHER an input samplesheet or a directory! Just pick one.' 
        } else { // if only samplesheet is passed check to make sure input is an actual file
            def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = null //keep input directory null if not passed
            // get full path for input and make channel
            if (params.input) { ch_input = file(params.input) }
        }
    } else {
        if (params.indir != null ) { // if no samplesheet is passed, but an input directory is given
            ch_input = null //keep samplesheet input null if not passed
            def checkPathParamList = [ params.indir, params.multiqc_config, params.kraken2db ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = Channel.fromPath(params.indir, relative: true, type: 'dir')
        } else { // if no samplesheet is passed and no input directory is given
            exit 1, 'For -entry UPDATE_CDC_PHOENIX: You need EITHER an input samplesheet or a directory!' 
        }
    }

    main:
        UPDATE_PHOENIX_WF ( ch_input, ch_input_indir, ch_versions )

    emit:
        mlst             = UPDATE_PHOENIX_WF.out.mlst
        amrfinder_output = UPDATE_PHOENIX_WF.out.amrfinder_output
        gamma_ar         = UPDATE_PHOENIX_WF.out.gamma_ar
        phx_summary      = UPDATE_PHOENIX_WF.out.phx_summary
        //output for phylophoenix
        griphin_tsv      = UPDATE_PHOENIX_WF.out.griphin_tsv
        griphin_excel    = UPDATE_PHOENIX_WF.out.griphin_excel
}

//
// WORKFLOW: Entry point for combining multiple griphin files
//
workflow COMBINE_GRIPHINS {
    // Check mandatory parameters
    ch_versions = Channel.empty() // Used to collect the software versions
    // Check input path parameters to see if they exist
    if (params.input != null ) {  // if a samplesheet is passed
        //input_samplesheet_path = Channel.fromPath(params.input, relative: true)
        if (params.indir != null ) { //if samplesheet is passed and an input directory exit
            exit 1, 'For -entry COMBINE_GRIPHINS: --indir is not a valid parameter, please pass a samplesheet and with --input.' 
        } else { // if only samplesheet is passed check to make sure input is an actual file
            def checkPathParamList = [ params.input, params.multiqc_config ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            // get full path for input and make channel
            if (params.input) { ch_input = file(params.input) }
            //check outdir
            if (params.outdir == "${launchDir}/phx_output" ) { 
                println("${orange}Warning: No outdir was passed, combined griphin file will be saved to the default ${launchDir}/phx_output.${reset}")
            } else {
                // Allow outdir to be relative
                outdir = Channel.fromPath(params.outdir, relative: true)
            }
        }
    } else {
        exit 1, 'For -entry COMBINE_GRIPHINS: --indir is not a valid parameter, please pass a samplesheet and with --input.' 
    }

    //no griphins to start - they should be in the input samplesheet
    input_griphins_ch = null
    //input_griphins_tsv_ch = null

    main:
        COMBINE_GRIPHINS_WF ( input_griphins_ch, ch_input, outdir, ch_versions )

}

/*
========================================================================================
    RUN Species specific WORKFLOWS - waiting for completed validation to be released with v2.3.0
========================================================================================
*/

//
// WORKFLOW: Entry point for running C. diff specific pipeline as standalone
//
workflow CENTAR {
    // comment out to run CENTAR 
    exit 1, "Sorry, -entry CENTAR hasn't completed its validation yet and will be released in another version of PHoeNIx!"

    // Check mandatory parameters
    ch_versions = Channel.empty() // Used to collect the software versions
    // Check input path parameters to see if they exist
    if (params.input != null ) {  // if a samplesheet is passed
        if (params.indir != null ) { //if samplesheet is passed and an input directory exit
            exit 1, 'For -entry RUN_CENTAR: You need EITHER an input samplesheet or a directory! Just pick one.' 
        } else { // if only samplesheet is passed check to make sure input is an actual file
            def checkPathParamList = [ params.input, params.multiqc_config ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            ch_input_indir = null //keep input directory null if not passed
            // get full path for input and make channel
            if (params.input) { ch_input = file(params.input) }
            // Allow outdir to be relative !!!! Does this need to be changed if outdir is empty??
            outdir = Channel.fromPath(params.outdir, relative: true)
            //griph_out = Channel.fromPath(params.griphin_out, relative: true)
        }
    } else {
        if (params.indir != null ) { // if no samplesheet is passed, but an input directory is given
            ch_input = null //keep samplesheet input null if not passed
            def checkPathParamList = [ params.indir, params.multiqc_config ]
            for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
            //make sure a directory is passed 
            if (new File(params.indir).isDirectory()){
                ch_input_indir = Channel.fromPath(params.indir, relative: true)
            } else {
                exit 1, 'You passed a file with --indir and a directory is required. Or use --input'
            }
            if (params.outdir == "${launchDir}/phx_output" ) { 
                outdir = params.indir
                println("${orange}Warning: No outdir was passed, so CENTAR files will be saved to the indir ${outdir}.${reset}")
            } else {
                // Allow outdir to be relative
                outdir = Channel.fromPath(params.outdir, relative: true)
                //griph_out = Channel.fromPath(params.griphin_out, relative: true)
            }
        } else { // if no samplesheet is passed and no input directory is given
            exit 1, 'For -entry CENTAR: You need EITHER an input samplesheet or a directory!' 
        }
    }

    //make sure outdir and griphin_out aren't passed at the same time
    if (params.griphin_out != "${launchDir}" && params.outdir != "${launchDir}/phx_output"){
        exit 1, "When using --outdir with CENTAR you can't use --griphin_out as --outdir directs all CENTAR and GRiPHin summary files to outdir. Please rerun with only one of these parameters." 
    }
    // check if the wgmlst_container was passed
    if (params.wgmlst_container == null) { println("${orange}Warning: No path was passed for --wgmlst_container so ribotyping will not be reported.${reset}") }

    main:
        RUN_CENTAR ( ch_input, ch_input_indir, ch_versions, outdir )

    emit:
        //output for phylophoenix
        griphins_excel   = RUN_CENTAR.out.griphins_excel
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {
    if(params.mode =="PHOENIX") {
        PHOENIX()
    } else if(params.mode =="CDC_PHOENIX") {
        CDC_PHOENIX()
    } else if(params.mode =="SRA") {
        SRA()
    } else if(params.mode =="CDC_SRA") {
        CDC_SRA()
    } else if(params.mode =="SCAFFOLDS") {
        SCAFFOLDS()
    } else if(params.mode =="CDC_SCAFFOLDS") {
        CDC_SCAFFOLDS()
    } else if(params.mode =="UPDATE_PHOENIX") {
        UPDATE_PHOENIX()
    } else if(params.mode =="COMBINE_GRIPHINS") {
        COMBINE_GRIPHINS()
    //} else if(params.mode =="CENTAR") {
    //    CENTAR()
    } else {
        exit 1, 'Please select an entry point either: PHOENIX, CDC_PHOENIX, SCAFFOLDS, CDC_SCAFFOLDS, SRA, CDC_SRA, UPDATE_PHOENIX and COMBINE_GRIPHINS'
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/