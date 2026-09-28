//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { SAMPLESHEET_CHECK as SAMPLESHEET_CHECK_2 } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    update_griphin

    main:
        if (update_griphin == false && params.mode_upper != "PHOENIX_LR" && params.mode_upper != "PHOENIX_HYBRID") {
           reads = SAMPLESHEET_CHECK ( samplesheet, true, false, false, false, [] ) // last [] used for --pipeline update_phoenix to get meta.full_project_id - to make sure things are published to the right dir in --input
                .csv.splitCsv ( header:true, sep:',' )
                .map { create_fastq_channels(it) }
            griphins = Channel.empty()
            long_read = [] //just an empty channel to keep it runnin
        } else if (params.mode_upper == "PHOENIX_LR") {
            SAMPLESHEET_CHECK ( samplesheet, false, false, false, false, [] )
                .csv
                .splitCsv ( header:true, sep:',' )
                .map { create_LR_fastq_channel(it) }
                .set { reads }
                griphins = Channel.empty()
                long_read = [] //just an empty channel to keep it runnin
        } else if (params.mode_upper == "PHOENIX_HYBRID") {
            SAMPLESHEET_CHECK ( samplesheet, false, false, false, false, [] )
                .csv
                .splitCsv ( header:true, sep:',' )
                .map { create_hybrid_fastq_channel(it) }
                .set { all_reads }

            // Split `all_reads` into two separate channels
            reads = all_reads.map { meta, reads -> [ [id:meta.id, single_end: false], reads[0..1] ] } // Take only fastq_1 and fastq_2, need single_end to be false for fastp
            long_read = all_reads.map { meta, reads -> [ meta, reads[2] ] } // Take only long_read
            griphins = Channel.empty()
        } else {
            griphins = SAMPLESHEET_CHECK ( samplesheet, false, false, false, true, [] ).csv.splitCsv(header:false, sep:',')
                .map{ row ->
                    if (!file(row[0]).exists()) { 
                        exit 1, "ERROR: Please check input samplesheet -> ${row[0]} does not exist!\n"
                    } else { return row }}.collect()
                .map{ rows ->
                    if (rows.size() < 2) {
                        exit 1, "ERROR: Need at least 2 rows in the griphin samplesheet, but found ${rows.size()}.\n"
                    } else { return rows }}
            reads = Channel.empty()
        }

    emit:
        griphins          = griphins
        reads             = reads                          // channel: [ val(meta), [ reads ] ]
        long_read         = long_read                      // channel: [ val(meta), [ long_read ] ]
        valid_samplesheet = SAMPLESHEET_CHECK.out.csv
        versions          = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

def check_griphins_exist(List row) {
    // Check if the file path exists
    if (!file(row[0]).exists()) { 
        exit 1, "ERROR: Please check input samplesheet -> ${row[0]} does not exist!\n"
    }
    return fileList
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

// Function to get list of [ meta, [ fastq_1, fastq_2, long_read ] ]
def create_hybrid_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (!file(row.fastq_2).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
    }
    if (!file(row.long_read).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> long read FastQ file does not exist!\n${row.long_read}"
    }
    array = [ meta, [ file(row.fastq_1), file(row.fastq_2), file(row.long_read) ] ]
    return array
}

// Function to get list of [ meta, [ scaffolds ] ]
def create_LR_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    if (!file(row.fastq).exists()) {
        exit 1, "ERROR: Please check input samplesheet ->  FastQ file does not exist!\n${row.fastq}"
    }
    return [ meta, file(row.fastq)]
}