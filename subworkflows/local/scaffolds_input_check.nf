//
// Check input samplesheet and get scaffolds channels
//

//create module below to address scaffolds
include { SCAFFOLDS_SAMPLESHEET_CHECK } from '../../modules/local/scaffolds_samplesheet_check'

workflow SCAFFOLDS_INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SCAFFOLDS_SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_scaff_channels(it) }
        .set { scaffolds }  //single fasta file not pairs

    emit:
    scaffolds                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions    // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_scaff_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    //meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.scaffolds_file).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Scaffolds file does not exist!\n${row.scaffolds_file}"
    } else {
    //if (meta.single_end) {
        array = [ meta, [ file(row.scaffolds_file) ] ]
    } //else {
        //if (!file(row.fastq_2).exists()) {
            //exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        //}
        //array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    //}
    return array
}
