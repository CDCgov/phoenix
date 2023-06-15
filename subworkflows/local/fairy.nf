//
// Check input samplesheet and get read channels
//

include { GET_RAW_STATS } from '../../modules/local/get_raw_stats'

workflow FAIRY_WKFLW {
    take:
    reads // file: /path/to/samplesheet.csv

    main: //RUN SCRIPT TO CHECK FILE INTEGRITY FOR SAMPLES
    /*SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channels(it) }
        .set { reads }*/
    //run GET_RAW_STATS
    //create new samplesheet w/only the passing FASTQs
    //create new fastq channel using the samplesheet check

    gzip -t ${reads[0]} 2>> ${fname}.txt
    gzip -t ${reads[1]} 2>> ${fnameB}.txt
    if grep -Fx "error" ${fname}.txt || grep -Fx "error" ${fnameB}.txt; then
        echo "FAIL" > ${prefix}_results.txt
    else
        echo "PASS" > ${prefix}_results.txt

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    valid_samplesheet = SAMPLESHEET_CHECK.out.csv
    versions          = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:] // ADD ADTL METADATA
    meta.id           = row.sample
    meta.single_end   = row.single_end.toBoolean()
    //meta.fairy        = //fairy

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