//
// workflow handles taking in either a samplesheet or directory and creates correct channels.
//

include { SAMPLESHEET_CHECK  } from '../../modules/local/samplesheet_check'
include { CREATE_SAMPLESHEET } from '../../modules/local/create_samplesheet'

workflow CREATE_SCAFFOLDS_INPUT_CHANNEL {
    take:
        indir        // params.indir
        samplesheet  // params.input
        ch_versions

    main:
        //if input directory is passed use it to gather assemblies otherwise use samplesheet
        if (indir != null) {
            if (params.scaffolds_ext !='.scaffolds.fa.gz') {
                def ext = ""
                if (!params.scaffolds_ext.toString().startsWith("*")) {
                    ext = "*" + params.scaffolds_ext.toString()
                } else {
                    ext = params.scaffolds_ext.toString()
                }
                def scaffolds_glob = append_to_path(params.indir.toString(), ext)
                //create scaffolds channel with meta information -- annoying, but you have to keep this in the brackets instead of having it once outside.
                scaffolds_ch = Channel.fromPath(scaffolds_glob) // use created regrex to get samples
                    .filter( it -> !(it =~ 'filtered') ) // remove samples that are *.filtered.scaffolds.fa.gz
                    .filter( it -> !(it =~ 'renamed') ) // remove samples that are *.renamed.scaffolds.fa.gz
                    .filter( it -> !(it =~ 'contig') ) // remove samples that are *.contigs.fa.gz
                    .map{ it -> create_meta(it, params.scaffolds_ext.toString())} // create meta for sample
                scaffolds_ch.view()
                // Checking regrex has correct extension
                scaffolds_ch.collect().map{ it -> check_scaffolds(it) }
            } else {
                def scaffolds_glob = append_to_path(params.indir.toString(),'*.scaffolds.fa.gz')
                //create scaffolds channel with meta information -- annoying, but you have to keep this in the brackets instead of having it once outside.
                scaffolds_ch = Channel.fromPath(scaffolds_glob) // use created regrex to get samples
                    .filter( it -> !(it =~ 'filtered') ) // remove samples that are *.filtered.scaffolds.fa.gz
                    .filter( it -> !(it =~ 'renamed') ) // remove samples that are *.renamed.scaffolds.fa.gz
                    .filter( it -> !(it =~ 'contig') ) // remove samples that are *.contigs.fa.gz
                    .map{ it -> create_meta(it, params.scaffolds_ext.toString())} // create meta for sample
                    //.ifEmpty(exit 1, "ERROR: Looks like there isn't assemblies in the folder you passed. PHoeNIx doesn't search recursively!\n") // this doesn't work for some reason. 
                // Checking regrex has correct extension
                scaffolds_ch.collect().map{ it -> check_scaffolds(it) }
            }

            //get valid samplesheet for griphin step in cdc_scaffolds
            CREATE_SAMPLESHEET (
                indir
            )
            ch_versions = ch_versions.mix(CREATE_SAMPLESHEET.out.versions)

            valid_samplesheet = CREATE_SAMPLESHEET.out.samplesheet
        } else if (samplesheet != null) {
            // if a samplesheet was passed then use that to create the channel
            scaffolds_ch = SAMPLESHEET_CHECK( samplesheet, false, true, false ) 
                .csv
                .splitCsv( header:true, sep:',' )
                .map { create_assembly_channel(it) }
            //get valid samplesheet for griphin step in cdc_scaffolds
            valid_samplesheet = SAMPLESHEET_CHECK.out.csv
        } else {
            exit 1, 'You need EITHER an input samplesheet or a directory!' 
        }

    emit:
        scaffolds_ch      = scaffolds_ch       // channel: [ meta, [ scaffolds_file ] ]
        valid_samplesheet = valid_samplesheet
        versions          = ch_versions

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

def create_meta(sample, file_extension){
    '''Creating meta: [[id:sample1, single_end:true], $PATH/sample1.scaffolds.fa.gz]'''
    sample_name_minus_path = sample.toString().split('/')[-1] // get the last string after the last backslash
    sample_name = sample_name_minus_path.replaceAll(file_extension, "") // remove file extention to get only sample name 
    def meta = [:] // create meta array
    meta.id = sample_name
    meta.single_end = 'true'
    array = [ meta, sample ]  //file() portion provides full path
    return array
}

// Function to get list of [ meta, [ scaffolds_file ] ]
def create_assembly_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = true.toBoolean()

    // add path(s) of the assembly file(s) to the meta map
    def assembly_meta = []
    if (!file(row.assembly).exists()) {
        exit 1, "ERROR: Please check assembly samplesheet -> Assembly scaffolds file does not exist!\n${row.assembly}"
    }
    assembly_meta = [ meta, file(row.assembly) ]
    return assembly_meta
}

def check_scaffolds(scaffold_channel) {
    if (scaffold_channel[1].toString().endsWith(".fasta.gz") or scaffold_channel[1].toString().endsWith(".fa.gz") ) {
        //If there is the correct ending just move along
    } else {
         exit 1, "ERROR: No scaffolds found. Either your scaffold regrex is off (scaffolds files should end in either '.fa.gz' or ''.fasta.gz') or the directory provided doesn't contain assembly files." 
    }
}