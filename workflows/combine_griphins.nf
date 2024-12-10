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

/*
========================================================================================
    IMPORT LOCAL MODULES
========================================================================================
*/

include { UPDATE_GRIPHIN                       } from '../modules/local/updater/update_griphin'
include { UPDATE_GRIPHIN as UPDATE_CDC_GRIPHIN } from '../modules/local/updater/update_griphin'

/*
========================================================================================
    IMPORT LOCAL SUBWORKFLOWS
========================================================================================
*/

include { CREATE_INPUT_CHANNELS          } from '../subworkflows/local/create_input_channels'

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS as COMBINE_CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS as CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS  } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main_centar'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/

def add_meta(input_ch) {
    def meta = [:] // create meta array
    meta.project_id = input_ch.getName().replaceAll("_GRiPHin.xlsx", "") // get file name without extention
    def array = [ meta, input_ch ]
    return array
}

def findGRiPHinSummaryFiles(path, extension) {
        def files = file(path[0]).listFiles() // List all files in the directory
        files.findAll { it.name.endsWith('_GRiPHin_Summary.'+ extension) } // Find all matching files
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow COMBINE_GRIPHINS_WF {
    take:
        input_griphins_excel_ch
        //input_griphins_tsv_ch
        ch_input
        outdir_path
        ch_versions

    main:
        if (input_griphins_excel_ch != null) { // for running at the end of another entry 

            // bring all the griphins into one channel and pass one at a time to the UPDATE_GRIPHIN process
            //update_griphin_ch = input_griphins_excel_ch.collect().combine(input_griphins_tsv_ch.collect()).combine(outdir_path)

            // combine griphin files, the new one just created and the old one that was found in the project dir. 
            UPDATE_GRIPHIN (
                input_griphins_excel_ch.collect(),
                //input_griphins_tsv_ch.collect(),
                outdir_path,
                outdir_path.toString().split('/')[-1], 
                //outdir_path.map{ dir -> dir.toString().split('/')[-1].replace("]","")},
                params.coverage
            )
            ch_versions = ch_versions.mix(UPDATE_GRIPHIN.out.versions)

            /// need to figure out how to do this on a per directory 
            // Collecting the software versions
            CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS (
                ch_versions.map{meta_file, version -> version},
                ch_versions.map{meta_file, version -> meta_file},
                ch_versions.map{meta_file, version -> meta_file.toString()}
            )

        } else { // for running as its own entry point

            //based on the location of the excel files get the tsv files.
            excel_griphins = Channel.from(ch_input).splitCsv().map{it -> findGRiPHinSummaryFiles(it, "xlsx")}
                .flatten()  // Flatten the list of lists into a single list of files
                .filter { it != null && it.exists() }.collect() // Ensure the files exist and are not null

            /*tsv_griphins = excel_griphins.flatten().map{it -> findGRiPHinSummaryFiles(it,"tsv")}
                .flatten()  // Flatten the list of lists into a single list of files
                .filter { it != null && it.exists() }.collect() // Ensure the files exist and are not null

            //tsv_griphins.view()*/

            // combine griphin files, the new one just created and the old one that was found in the project dir. 
            UPDATE_GRIPHIN (
                excel_griphins, 
                outdir_path, 
                outdir_path.map{ dir -> dir.toString().split('/')[-1].replace("]","")},
                params.coverage
            )
            ch_versions = ch_versions.mix(UPDATE_GRIPHIN.out.versions)
            
            /// need to figure out how to do this on a per directory 
            // Collecting the software versions
            COMBINE_CUSTOM_DUMPSOFTWAREVERSIONS (
                ch_versions.unique().collectFile(name: 'collated_versions.yml')
            )
        }

    emit:
        //output for phylophoenix
        griphin_tsv      = UPDATE_GRIPHIN.out.griphin_tsv_report
        griphin_excel    = UPDATE_GRIPHIN.out.griphin_report
        //dir_samplesheet  = UPDATE_CDC_GRIPHIN.out.converted_samplesheet*/
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

// Adding if/else for running on ICA
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
*/
