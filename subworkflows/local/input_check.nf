//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK    } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv
        entry_point // if LR we need to return a single fastq, if Illumina we need to return two fastq files

    main:
        if (entry_point == "Illumina") {
            SAMPLESHEET_CHECK ( samplesheet, true, false, false, false )
                .csv
                .splitCsv ( header:true, sep:',' )
                .map { create_fastq_channels(it) }
                .set { reads }
        } else if (entry_point == "Nanopore") {
            SAMPLESHEET_CHECK ( samplesheet, false, false, false, true )
                .csv
                .splitCsv ( header:true, sep:',' )
                .map { create_LR_fastq_channel(it) }
                .set { reads }
        } else{
            exit 1, "ERROR: entry_point variable needs to be set to Nanopore or Illumina\n"
        }

    emit:
        reads                                              // channel: [ val(meta), [ reads ] ]
        valid_samplesheet = SAMPLESHEET_CHECK.out.csv
        versions          = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}

// Function to get list of [ meta, [ scaffolds ] ]
def create_LR_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.sample
    meta.single_end = row.single_end

    // add path(s) of the fastq file(s) to the meta map
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet ->  FastQ file does not exist!\n${row.fastq}"
    }
    return [ meta, file(row.fastq)]
}
