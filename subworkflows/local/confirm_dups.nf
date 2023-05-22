
//
// workflow handles taking in either a samplesheet or directory and creates correct channels.
//

workflow CONFIRM_DUPS {
    take:
        combined_sra_ch
        metadata_csvs

    main:

    // Check for duplicate names
    confirmed_duplicates = metadata_csvs.map{meta, metadata_csv -> [ metadata_csv ] }.collect().map{metadata_csvs -> check_for_dups( metadata_csvs ) }

    // Use result of duplicate check to create new channel
    rename_sra_ch = combined_sra_ch.combine(confirmed_duplicates).map{meta, reads, metadata_csv, dups -> create_deduped_ch(meta, reads, metadata_csv, dups) }

    emit:
        rename_sra_out_ch = rename_sra_ch 
}

def create_deduped_ch(old_meta, reads, metadata_csv, dups){
    if (dups == "no") { // There are no duplicate sample names so use those for the meta
        // Create tuple with sample name in it for the tag in RENAME_SRA_FASTA
        // function is used to swap out SRR meta.id for the sample name from SRA
        def meta = [:] // create meta array
        meta.id = metadata_csv.readLines().get(1).split(',')[29].replaceAll(" ", "_").replaceAll("/", "_") // This gives the metadata sample name from the SRA, also some clean up
        output_array = [ meta, reads]
    } else {// There are duplicate sample names so use SRR number for naming
        // Just drop metadata_csv from the channel
        output_array = [ old_meta, reads]
    }
    return output_array
}

def check_for_dups(metadata_csvs) {
    def list_of_sample_names = []
    metadata_csvs.each {  // loop through each metadata_csv in the list
        sample_name = it.readLines().get(1).split(',')[29].replaceAll(" ", "_").replaceAll("/", "_") // This gives the metadata sample name from the SRA, also some clean up
        list_of_sample_names.add(sample_name) // add each sample name to list
    }
    // get length of original list of names
    int list_len = list_of_sample_names.size
    //get only unique names
    int unique_len = list_of_sample_names.unique().size
    def duplicates = unique_len < list_len? "yes": "no"
    return duplicates
}