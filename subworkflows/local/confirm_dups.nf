
//
// workflow handles taking in either a samplesheet or directory and creates correct channels.
//

workflow CONFIRM_DUPS {
    take:
        combined_sra_ch
        metadata_csvs

    main:

    if (params.use_sra) { // if --use_srr is passed then this overrides attempts to have samples named by sample names
        rename_sra_ch = combined_sra_ch.map{meta, reads, metadata_csv-> [ meta, reads] }
    } else {
        // Check for duplicate names
        confirmed_duplicates = metadata_csvs.map{meta, metadata_csv -> [ metadata_csv ] }.collect().map{metadata_csvs -> check_for_dups( metadata_csvs ) }
        // Use result of duplicate check to create new channel
        rename_sra_ch = combined_sra_ch.combine(confirmed_duplicates).map{meta, reads, metadata_csv, dups -> create_deduped_ch(meta, reads, metadata_csv, dups) }
    }

    emit:
        rename_sra_out_ch = rename_sra_ch 
}

def create_deduped_ch(old_meta, reads, metadata_csv, dups){
    if (dups == "no") { // There are no duplicate sample names so use those for the meta
        // Create tuple with sample name in it for the tag in RENAME_SRA_FASTA
        // function is used to swap out SRR meta.id for the sample name from SRA
        def meta = [:] // create meta array
        try {
            def lines = metadata_csv.readLines()
            // Sanity checks
            if (lines.size() > 1) {
                def cols = lines[1].split(',')
                if (cols.size() > 29 && cols[29]) {
                    meta.id = cols[29].replaceAll(" ", "_").replaceAll("/", "_")
                } else {
                    throw new IndexOutOfBoundsException("Missing SampleName column")
                }
            } else {
                throw new IndexOutOfBoundsException("File has no data rows")
            }
        } catch (Exception e) {
            // Fallback to filename-based sample name
            meta.id = metadata_csv.toString().replace("_sra_metadata.csv", "").tokenize("/")[-1]
        }
        //meta.id = metadata_csv.readLines().get(1).split(',')[29].replaceAll(" ", "_").replaceAll("/", "_") // This gives the metadata sample name from the SRA, also some clean up
        output_array = [ meta, reads]
    } else {// There are duplicate sample names so use SRR number for naming
        // Just drop metadata_csv from the channel
        output_array = [ old_meta, reads]
    }
    return output_array
}

def check_for_dups(metadata_csvs) {
    def list_of_sample_names = []
    metadata_csvs.each { file -> // loop through each metadata_csv in the list
        def sample_name
        try {
            def lines = file.readLines()
            // Sanity checks
            if (lines.size() > 1) {
                def cols = lines[1].split(',')
                if (cols.size() > 29 && cols[29]) {
                    sample_name = cols[29].replaceAll(" ", "_").replaceAll("/", "_")
                } else {
                    throw new IndexOutOfBoundsException("Missing SampleName column")
                }
            } else {
                throw new IndexOutOfBoundsException("File has no data rows")
            }
        } catch (Exception e) {
            // Fallback to filename-based sample name
            sample_name = file.toString().replace("_sra_metadata.csv", "").tokenize("/")[-1]
        }
        list_of_sample_names << sample_name
    }
    def list_len = list_of_sample_names.size()
    def unique_len = list_of_sample_names.unique().size()
    return unique_len < list_len ? "yes" : "no"
}
/* Old version of check_for_dups function for reference
        sample_name = it.readLines().get(1).split(',')[29].replaceAll(" ", "_").replaceAll("/", "_") // This gives the metadata sample name from the SRA, also some clean up
        sample_name = it.toString().split("/")[-1].replace("_sra_metadata.csv","")
        list_of_sample_names.add(sample_name) // add each sample name to list
    }
    // get length of original list of names
    int list_len = list_of_sample_names.size
    //get only unique names
    int unique_len = list_of_sample_names.unique().size
    def duplicates = unique_len < list_len? "yes": "no"
    return duplicates
}*/