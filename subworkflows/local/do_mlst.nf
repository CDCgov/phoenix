//
// Subworkflow: Running srst2_MLST
//

include { MLST                           } from '../../modules/local/mlst'
include { SRST2_MLST                     } from '../../modules/local/srst2_mlst'
include { GET_MLST_SRST2                 } from '../../modules/local/get_mlst_srst2'
include { CHECK_MLST                     } from '../../modules/local/check_mlst'
include { CHECK_MLST as CHECK_ALL_MLST   } from '../../modules/local/check_mlst_external'

workflow DO_MLST {
    take:
        trimmed_assembly  // channel: tuple val(meta), path(assembly)
        paired_reads      // channel: tuple val(meta), path(reads), path(paired_reads):FASTP_TRIMD.out.reads.map
        taxonomy          // channel: tuple val(meta), path(taxonomy): DETERMINE_TAXA_ID.out.taxonomy
        do_srst2_mlst     //

    main:
        ch_versions = Channel.empty() // Used to collect the software versions

        mlst_ch = trimmed_assembly.map{meta, fasta    -> [[id:meta.id], fasta]}

        MLST (
            mlst_ch
        )
        ch_versions = ch_versions.mix(MLST.out.versions)

        check_main_mlst_ch = MLST.out.tsv.map{     meta, tsv      -> [[id:meta.id], tsv]}\
        .join(taxonomy.map{meta, taxonomy -> [[id:meta.id], taxonomy]}, by: [0])\

        // Combining and adding flare to all MLST outputs
        CHECK_MLST (
            check_main_mlst_ch
        )
        ch_versions = ch_versions.mix(CHECK_MLST.out.versions)




        // if cdc_phoenix, we want to do srst2, if it needs it, otherwise just skip
        if (do_srst2_mlst == true) {
            pre_GET_MLST_SRST2_ch = taxonomy.map{meta, taxonomy    -> [[id:meta.id], taxonomy]}\
            .join(CHECK_MLST.out.status.splitCsv(strip:true).map{   meta, status -> [[id:meta.id], status]},  by: [0])

            // Runs the getMLST portion of the srst2 mlst script to find right scheme to compare against
            GET_MLST_SRST2 (
                pre_GET_MLST_SRST2_ch
            )
            ch_versions = ch_versions.mix(GET_MLST_SRST2.out.versions)

            // Combining weighted kraken report with the FastANI hit based on meta.id
            mid_srst2_ch = paired_reads.map{meta, reads    -> [[id:meta.id], reads]}\
            .join(GET_MLST_SRST2.out.getMLSTs.map{   meta, getMLSTs -> [[id:meta.id], getMLSTs]},  by: [0])\
            .join(GET_MLST_SRST2.out.fastas.map{     meta, fastas   -> [[id:meta.id], fastas]},    by: [0])\
            .join(GET_MLST_SRST2.out.profiles.map{   meta, profiles -> [[id:meta.id], profiles]},  by: [0])\
            .join(CHECK_MLST.out.status.splitCsv(strip:true).map{   meta, status -> [[id:meta.id], status]},  by: [0])

            // Identifying mlst genes in trimmed reads
            SRST2_MLST (
                mid_srst2_ch
            )
            ch_versions = ch_versions.mix(SRST2_MLST.out.versions)

            combined_mlst_ch = MLST.out.tsv.map{     meta, tsv           -> [[id:meta.id], tsv]}\
            .join(SRST2_MLST.out.mlst_results_temp.map{   meta, mlst_results  -> [[id:meta.id], mlst_results]}, by: [0])\
            .join(taxonomy.map{meta, taxonomy      -> [[id:meta.id], taxonomy]},     by: [0])\
            .join(SRST2_MLST.out.empty_checker.splitCsv(strip:true).map{   meta, status -> [[id:meta.id], status]},  by: [0])
                
            // Combining and adding flare to all MLST outputs
            CHECK_ALL_MLST (
                combined_mlst_ch
            )
            ch_versions = ch_versions.mix(CHECK_ALL_MLST.out.versions)

            checked_mlst_ch = CHECK_ALL_MLST.out.checked_MLSTs.map{meta, checked_MLSTs -> [ [id:meta.id], checked_MLSTs]}

            // Defining out channel
            //raw_check_channel = CHECK_ALL_MLST.out.checked_MLSTs.ifempty(CHECK_MLST.out.checked_MLSTs)
            //checked_mlst_ch = raw_check_channel.map{meta, checked_MLSTs -> [ [id:meta.id], checked_MLSTs]}.ifEmpty(CHECK_MLST.out.checked_MLSTs.map{meta, checked_MLSTs -> [ [id:meta.id], checked_MLSTs]})
            checked_mlst_ch = CHECK_ALL_MLST.out.checked_MLSTs.map{meta, checked_MLSTs -> [ [id:meta.id], checked_MLSTs]}
        } else {
            checked_mlst_ch = CHECK_MLST.out.checked_MLSTs.map{meta, checked_MLSTs -> [ [id:meta.id], checked_MLSTs]}
        }

    emit:
        checked_mlsts               = checked_mlst_ch
        versions                    = ch_versions // channel: [ versions.yml ]
}

