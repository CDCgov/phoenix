//
// Subworkflow: run Kraken2
//

include { KRAKEN2_KRAKEN2                               } from '../../modules/local/kraken2'
include { KRAKEN2_KRONA                                 } from '../../modules/local/krakentools_kreport2krona'
include { KRONA_KTIMPORTTEXT                            } from '../../modules/local/ktimporttext'
include { KRAKENTOOLS_KREPORT2MPA                       } from '../../modules/local/krakentools_kreport2mpa'
include { KRAKENTOOLS_MAKEKREPORT                       } from '../../modules/local/krakentools_makekreport'
include { KRAKEN_BEST_HIT                               } from '../../modules/local/kraken_bh'

workflow KRAKEN2_WF {
    take:
    fasta           // channel: tuple (meta) path(read_R1, reads_R2) or tuple (meta) path(scaffolds)
    fairy_check     // GET_RAW_STATS.out.outcome or SCAFFOLD_COUNT_CHECK.out.outcome
    type            // val: trimd, asmbld or wtasmbld 
    qc_stats        //GATHERING_READ_QC_STATS.out.fastp_total_qc
    quast           //QUAST.out.report_tsv --> only for wtasmbld and asmbld
    kraken2_db_path
    seq_type        // either "scaffolds" or "reads"

    main:
    ch_versions     = Channel.empty() // Used to collect the software versions
    // Add in krakendb into the fasta channel so each fasta has a krakendb to go with it. If you don't do this then only one sample goes through pipeline
    if (type=="trimd") {
        // add in fairy to confirm reads are uncorrupted and correct
        fasta_ch = fasta.join(fairy_check.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [meta, [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0,0])\
        .combine(kraken2_db_path)
    } else if(type=="asmbld" || type=="wtasmbld") {
        // add in scaffold_count_check so its confirmed that there are scaffolds in file post filtering
        if (seq_type =="scaffolds") {
            fasta_ch = fasta.join(fairy_check.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [[id:meta.id, single_end:meta.single_end], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0])\
            .combine(kraken2_db_path)
        } else if(seq_type =="reads"){
            fasta_ch = fasta.join(fairy_check.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [[id:meta.id, single_end:true], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0])\
            .combine(kraken2_db_path)
        }
    }

    // Checking for Contamination in trimmed reads
    KRAKEN2_KRAKEN2 (
        fasta_ch, type, params.save_output_fastqs, params.save_reads_assignment
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

    // Create mpa file
    KRAKENTOOLS_KREPORT2MPA (
        KRAKEN2_KRAKEN2.out.report
    )
    ch_versions = ch_versions.mix(KRAKENTOOLS_KREPORT2MPA.out.versions)

    if (type == "trimd" || type == "asmbld"){
        // Converting kraken report to krona file to have hierarchical output in krona plot
        KRAKEN2_KRONA (
            KRAKEN2_KRAKEN2.out.report, type
        )
        report = KRAKEN2_KRAKEN2.out.report
    } else if (type == "wtasmbld"){        
        // Add in krakendb into the kraken reads channel so each fasta has a krakendb to go with it. 
        make_report_ch = KRAKEN2_KRAKEN2.out.classified_reads_assignment.combine(kraken2_db_path)

        // Create weighted kraken report based on scaffold length
        KRAKENTOOLS_MAKEKREPORT (
            make_report_ch
        )
        ch_versions = ch_versions.mix(KRAKENTOOLS_MAKEKREPORT.out.versions)

        // Converting kraken report to krona file to have hierarchical output in krona plot
        KRAKEN2_KRONA (
            KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report, "wtasmbld"
        )
        
        report = KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report
    }
    ch_versions = ch_versions.mix(KRAKEN2_KRONA.out.versions)

    
    // Create krona plot from kraken report
    KRONA_KTIMPORTTEXT (
        KRAKEN2_KRONA.out.krona, type
    )
    ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT.out.versions)

    if (type == "trimd"){
        // Combining kraken report with quast report based on meta.id
        kraken_bh_ch = KRAKEN2_KRAKEN2.out.report.map{meta, report         -> [[id:meta.id], report]}\
        .join(qc_stats.map{ meta, fastp_total_qc -> [[id:meta.id], fastp_total_qc]}, by: [0])
    } else if (type == "asmbld"){
        // Combining kraken report with quast report based on meta.id
        kraken_bh_ch = KRAKEN2_KRAKEN2.out.report.map{meta, report     -> [[id:meta.id], report]}\
        .join(quast.map{ meta, report_tsv -> [[id:meta.id], report_tsv]}, by: [0])
    } else if (type == "wtasmbld"){
        // Combining kraken report with quast report based on meta.id
        kraken_bh_ch = KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report.map{meta, kraken_weighted_report -> [[id:meta.id], kraken_weighted_report]}\
        .join(quast.map{                                                               meta, report_tsv             -> [[id:meta.id], report_tsv]}, by: [0])
    }
        
    // Getting Kraken best hit for assembled data
    KRAKEN_BEST_HIT (
        kraken_bh_ch, type
    )
    ch_versions = ch_versions.mix(KRAKEN_BEST_HIT.out.versions)

    emit:
    report        = report
    k2_bh_summary = KRAKEN_BEST_HIT.out.ksummary
    krona_html    = KRONA_KTIMPORTTEXT.out.html
    versions      = ch_versions // channel: [ versions.yml ]
}