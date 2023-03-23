//
// workflow handles taking in either a samplesheet or directory and creates correct channels.
//

include { SCAFFOLDS_SAMPLESHEET_CHECK } from '../../modules/local/scaffolds_samplesheet_check'
include { CREATE_SAMPLESHEET          } from '../../modules/local/create_samplesheet'

workflow CREATE_INPUT_CHANNEL {
    take:
        indir        // params.indir
        samplesheet  // params.input

    main:
        //if input directory is passed use it to gather assemblies otherwise use samplesheet
        if (indir != null) {
            if (params.scaffold_regrex !='*/*/*.scaffolds.fa.gz') {
                def scaffolds_glob = append_to_path(params.indir.toString(), params.scaffold_regrex)
                //create scaffolds channel with meta information -- annoying, but you have to keep this in the brackets instead of having it once outside.
                scaffolds_ch = Channel.fromPath(scaffolds_glob) // use created regrex to get samples
                    .filter( it -> !(it =~ 'filtered') ) // remove samples that are *.filtered.scaffolds.fa.gz
                    .filter( it -> !(it =~ 'renamed') ) // remove samples that are *.renamed.scaffolds.fa.gz
                    .map{ it -> create_meta(it)} // create meta for sample
            } else {
                def scaffolds_glob = append_to_path(params.indir.toString(),'*/*/*.scaffolds.fa.gz')
                //create scaffolds channel with meta information -- annoying, but you have to keep this in the brackets instead of having it once outside.
                scaffolds_ch = Channel.fromPath(scaffolds_glob) // use created regrex to get samples
                    .filter( it -> !(it =~ 'filtered') ) // remove samples that are *.filtered.scaffolds.fa.gz
                    .filter( it -> !(it =~ 'renamed') ) // remove samples that are *.renamed.scaffolds.fa.gz
                    .map{ it -> create_meta(it)} // create meta for sample
            }

            //get valid samplesheet for griphin step in cdc_scaffolds
            CREATE_SAMPLESHEET (
                indir
            )

            valid_samplesheet = CREATE_SAMPLESHEET.out.samplesheet
        } else if (samplesheet != null) {
            // if a samplesheet was passed then use that to create the channel
            scaffolds_ch = SCAFFOLDS_SAMPLESHEET_CHECK( samplesheet ) 
                .csv
                .splitCsv( header:true, sep:',' )
                .map { create_assembly_channel(it) }
            //get valid samplesheet for griphin step in cdc_scaffolds
            valid_samplesheet = SCAFFOLDS_SAMPLESHEET_CHECK.out.csv
        } else {
            exit 1, 'You need EITHER an input samplesheet or a directory!' 
        }

    emit:
        scaffolds_ch      = scaffolds_ch       // channel: [ meta, [ scaffolds_file ] ]
        valid_samplesheet = valid_samplesheet

}

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/

def append_to_path(full_path, string) {
    if (full_path.toString().endsWith('/')) {
        new_string = full_path.toString() + string
    }  else {
        new_string = full_path.toString() + '/' + string
    }
    return new_string
}

def create_meta(sample){
    '''Creating meta: [[id:sample1, single_end:true], $PATH/sample1.scaffolds.fa.gz]'''
    full_sample_name = sample.toString().split('/')[-1] // get the last string after the last backslash
    sample_name = full_sample_name.replaceAll(".scaffolds.fa.gz", "") // remove _samplesheet.csv from path 
    def meta = [:] // create meta array
    meta.id = sample_name
    meta.single_end = 'true'
    array = [ meta, file(sample) ]  //file() portion provides full path
    return array
}

// Function to get list of [ meta, [ scaffolds_file ] ]
def create_assembly_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = 'true'

    // add path(s) of the assembly file(s) to the meta map
    def assembly_meta = []
    if (!file(row.assembly).exists()) {
        exit 1, "ERROR: Please check assembly samplesheet -> Assembly scaffolds file does not exist!\n${row.assembly}"
    }
    assembly_meta = [ meta, file(row.assembly) ]
    return assembly_meta
}