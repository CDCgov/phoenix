//
// Subworkflow: run Kraken2
//

include { KRAKEN2_KRAKEN2 as KRAKEN2_TRIMD                  } from '../../modules/local/kraken2'
include { KRAKEN2_KRAKEN2 as KRAKEN2_ASMBLD                 } from '../../modules/local/kraken2'
include { KRAKEN2_KRAKEN2 as KRAKEN2_WTASMBLD               } from '../../modules/local/kraken2'
include { KRAKEN2_KRONA as KREPORT2KRONA_TRIMD              } from '../../modules/local/krakentools_kreport2krona'
include { KRAKEN2_KRONA as KREPORT2KRONA_ASMBLD             } from '../../modules/local/krakentools_kreport2krona'
include { KRAKEN2_KRONA as KREPORT2KRONA_WTASMBLD           } from '../../modules/local/krakentools_kreport2krona'
include { KRONA_KTIMPORTTEXT as KRONA_KTIMPORTTEXT_TRIMD    } from '../../modules/local/ktimporttext'
include { KRONA_KTIMPORTTEXT as KRONA_KTIMPORTTEXT_ASMBLD   } from '../../modules/local/ktimporttext'
include { KRONA_KTIMPORTTEXT as KRONA_KTIMPORTTEXT_WTASMBLD } from '../../modules/local/ktimporttext'
include { KRAKENTOOLS_KREPORT2MPA as KREPORT2MPA_TRIMD      } from '../../modules/local/krakentools_kreport2mpa'
include { KRAKENTOOLS_KREPORT2MPA as KREPORT2MPA_ASMBLD     } from '../../modules/local/krakentools_kreport2mpa'
include { KRAKENTOOLS_MAKEKREPORT                           } from '../../modules/local/krakentools_makekreport'
include { KRAKEN_BEST_HIT as KRAKEN2_BH_TRIMD               } from '../../modules/local/kraken_bh'
include { KRAKEN_BEST_HIT as KRAKEN2_BH_ASMBLD              } from '../../modules/local/kraken_bh'
include { KRAKEN_BEST_HIT as KRAKEN2_BH_WTASMBLD            } from '../../modules/local/kraken_bh'

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

    if(type =="trimd") {

        // Checking for Contamination in trimmed reads
        KRAKEN2_TRIMD (
            fasta_ch, "trimd", true, true
        )
        ch_versions = ch_versions.mix(KRAKEN2_TRIMD.out.versions)

        // Create mpa file
        KREPORT2MPA_TRIMD (
            KRAKEN2_TRIMD.out.report
        )
        ch_versions = ch_versions.mix(KREPORT2MPA_TRIMD.out.versions)

        // Converting kraken report to krona file to have hierarchical output in krona plot
        KREPORT2KRONA_TRIMD (
            KRAKEN2_TRIMD.out.report, "trimd"
        )
        ch_versions = ch_versions.mix(KREPORT2KRONA_TRIMD.out.versions)

        // Create krona plot from kraken report
        KRONA_KTIMPORTTEXT_TRIMD (
            KREPORT2KRONA_TRIMD.out.krona, "trimd"
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT_TRIMD.out.versions)

        // Combining kraken report with quast report based on meta.id
        kraken_bh_trimd_ch = KRAKEN2_TRIMD.out.report.map{meta, report         -> [[id:meta.id], report]}\
        .join(qc_stats.map{                               meta, fastp_total_qc -> [[id:meta.id], fastp_total_qc]}, by: [0])

        // Getting Kraken best hit for assembled data
        KRAKEN2_BH_TRIMD (
            kraken_bh_trimd_ch, "trimd"
        )
        ch_versions = ch_versions.mix(KRAKEN2_BH_TRIMD.out.versions)

        report        = KRAKEN2_TRIMD.out.report
        k2_bh_summary = KRAKEN2_BH_TRIMD.out.ksummary
        krona_html    = KRONA_KTIMPORTTEXT_TRIMD.out.html

    } else if(type =="asmbld") {

        // Checking for Contamination in scaffolds
        KRAKEN2_ASMBLD (
            fasta_ch, "asmbld", true, true
        )
        ch_versions = ch_versions.mix(KRAKEN2_ASMBLD.out.versions)

        // Create mpa file
        KREPORT2MPA_ASMBLD (
            KRAKEN2_ASMBLD.out.report
        )
        ch_versions = ch_versions.mix(KREPORT2MPA_ASMBLD.out.versions)

        // Converting kraken report to krona file to have hierarchical output in krona plot
        KREPORT2KRONA_ASMBLD (
            KRAKEN2_ASMBLD.out.report, "asmbld"
        )
        ch_versions = ch_versions.mix(KREPORT2KRONA_ASMBLD.out.versions)

        // Create krona plot from kraken report
        KRONA_KTIMPORTTEXT_ASMBLD (
            KREPORT2KRONA_ASMBLD.out.krona, "asmbld"
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT_ASMBLD.out.versions)

        kraken_bh_asmbld_ch = KRAKEN2_ASMBLD.out.report.map{meta, report     -> [[id:meta.id], report]}\
        .join(quast.map{                                    meta, report_tsv -> [[id:meta.id], report_tsv]}, by: [0])

        // Getting Kraken best hit for assembled data
        KRAKEN2_BH_ASMBLD (
            kraken_bh_asmbld_ch, "asmbld"
        )
        ch_versions = ch_versions.mix(KRAKEN2_BH_ASMBLD.out.versions)

        report        = KRAKEN2_ASMBLD.out.report
        k2_bh_summary = KRAKEN2_BH_ASMBLD.out.ksummary
        krona_html    = KRONA_KTIMPORTTEXT_ASMBLD.out.html

    } else if(type=="wtasmbld") {

        // Getting species ID as back up for FastANI and checking contamination isn't in assembly
        KRAKEN2_WTASMBLD (
            fasta_ch, "wtasmbld", true, true
        )
        ch_versions = ch_versions.mix(KRAKEN2_WTASMBLD.out.versions)

        // Add in krakendb into the kraken reads channel so each fasta has a krakendb to go with it. 
        make_report_ch = KRAKEN2_WTASMBLD.out.classified_reads_assignment.combine(kraken2_db_path)

        // Create weighted kraken report based on scaffold length
        KRAKENTOOLS_MAKEKREPORT (
            make_report_ch
        )
        ch_versions = ch_versions.mix(KRAKENTOOLS_MAKEKREPORT.out.versions)

        // Converting kraken report to krona file to have hierarchical output in krona plot
        KREPORT2KRONA_WTASMBLD (
            KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report, "wtasmbld"
        )
        ch_versions = ch_versions.mix(KREPORT2KRONA_WTASMBLD.out.versions)

        // Combining kraken report with quast report based on meta.id
        kraken_bh_wtasmbld_ch = KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report.map{meta, kraken_weighted_report -> [[id:meta.id], kraken_weighted_report]}\
        .join(quast.map{                                                               meta, report_tsv             -> [[id:meta.id], report_tsv]}, by: [0])

        // Getting Kraken best hit for assembled data
        KRAKEN2_BH_WTASMBLD (
            kraken_bh_wtasmbld_ch, "wtasmbld"
        )
        ch_versions = ch_versions.mix(KRAKEN2_BH_WTASMBLD.out.versions)

        KRONA_KTIMPORTTEXT_WTASMBLD (
            KREPORT2KRONA_WTASMBLD.out.krona, "wtasmbld"
        )
        ch_versions = ch_versions.mix(KRONA_KTIMPORTTEXT_WTASMBLD.out.versions)

        report        = KRAKENTOOLS_MAKEKREPORT.out.kraken_weighted_report
        k2_bh_summary = KRAKEN2_BH_WTASMBLD.out.ksummary
        krona_html    = KRONA_KTIMPORTTEXT_WTASMBLD.out.html

    } else {
        println("Type options are: wtasmbld, asmbld or trimd")
    }

    emit:
    report        = report
    k2_bh_summary = k2_bh_summary
    krona_html    = krona_html
    versions      = ch_versions // channel: [ versions.yml ]

}