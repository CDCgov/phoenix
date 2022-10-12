#!/usr/bin/env nextflow
/*
========================================================================================
    CDCgov/phoenix
========================================================================================
    Github : https://github.com/CDCgov/phoenix
    Slack  : https://staph-b-dev.slack.com/channels/phoenix-h-dev
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { PHOENIX_EXTERNAL } from './workflows/phoenix'
include { PHOENIX_EXQC     } from './workflows/cdc_phoenix'
include { SRA_PHOENIX      } from './workflows/sra_phoenix'

//
// WORKFLOW: Run main cdcgov/phoenix analysis pipeline
//
workflow PHOENIX {
    PHOENIX_EXTERNAL ()
}

//
// WORKFLOW: Run internal version of cdcgov/phoenix analysis pipeline that includes BUSCO, SRST2 and KRAKEN_ASMBLED
//
workflow CDC_PHOENIX {
    PHOENIX_EXQC ()
}

/*
========================================================================================
    RUN SRA WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow SRA_LIST {
    SRA_PHOENIX ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
