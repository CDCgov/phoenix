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

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

//Check coverage is above its threshold
if (params.coverage.toInteger() < 30) { exit 1, 'The minimum coverage allowed for QA/QC purposes is 30 and is the default. Please choose a value >=30.' }
//Check path of kraken2db
if (params.kraken2db == null) { exit 1, 'Input path to kraken2db not specified!' }

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

//
// WORKFLOW: Run main cdcgov/phoenix analysis pipeline
//
workflow PHOENIX {
    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db] //removed , params.fasta to stop issue w/connecting to aws and igenomes not used
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters

    //input on command line
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'For -entry PHOENIX: Input samplesheet not specified!' }
    ch_versions = Channel.empty() // Used to collect the software versions

    main:
        PHOENIX_EXTERNAL ( ch_input, ch_versions, true )
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
    // Validate input parameters
    // Check input path parameters to see if they exist
    def checkPathParamList = [ params.input, params.multiqc_config, params.kraken2db]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    // Check mandatory parameters

    //input on command line
    if (params.input) { ch_input = file(params.input) } else { exit 1, 'For -entry CDC_PHOENIX: Input samplesheet not specified!' }
    ch_versions = Channel.empty() // Used to collect the software versions

    main:
        PHOENIX_EXQC ( ch_input, ch_versions, true )

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
        PHOENIX_EXTERNAL ( SRA_PREP.out.samplesheet, SRA_PREP.out.versions, false )

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
        PHOENIX_EXQC ( SRA_PREP.out.samplesheet, SRA_PREP.out.versions, false )

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
        SCAFFOLDS_EXTERNAL ( ch_input, ch_input_indir )

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
        SCAFFOLDS_EXQC ( ch_input, ch_input_indir )

    emit:
        scaffolds        = SCAFFOLDS_EXQC.out.scaffolds
        mlst             = SCAFFOLDS_EXQC.out.mlst
        amrfinder_output = SCAFFOLDS_EXQC.out.amrfinder_output
        gamma_ar         = SCAFFOLDS_EXQC.out.gamma_ar
        phx_summary      = SCAFFOLDS_EXQC.out.phx_summary
}

//
// WORKFLOW: Entry point for CLIA analysis
//

workflow CLIA {
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

    /*emit:
        scaffolds        = CLIA_INTERNAL.out.scaffolds
        trimmed_reads    = CLIA_INTERNAL.out.trimmed_reads
        amrfinder_report = CLIA_INTERNAL.out.amrfinder_report
        summary_report   = CLIA_INTERNAL.out.summary_report*/
}


//
// WORKFLOW: Entry point for updating phoenix mlst and ar output
//

workflow UPDATE_PHOENIX {

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
            exit 1, 'For -entry UPDATE_PHOENIX: You need EITHER an input samplesheet or a directory! Just pick one.' 
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
            exit 1, 'For -entry UPDATE_PHOENIX: You need EITHER an input samplesheet or a directory!' 
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

/*
========================================================================================
    THE END
========================================================================================
*/
