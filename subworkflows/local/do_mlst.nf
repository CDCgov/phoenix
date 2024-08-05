//
// Subworkflow: Running srst2_MLST
//

include { MLST                           } from '../../modules/local/mlst'
include { SRST2_MLST                     } from '../../modules/local/srst2_mlst_local'
include { GET_MLST_SRST2                 } from '../../modules/local/get_mlst_srst2'
include { CHECK_MLST                     } from '../../modules/local/check_mlst'
include { CHECK_MLST_WITH_SRST2          } from '../../modules/local/check_mlst_with_srst2'

workflow DO_MLST {
    take:
        trimmed_assembly     // channel: tuple val(meta), path(assembly): BBMAP_REFORMAT.out.filtered_scaffolds
        scaffold_count_check // SCAFFOLD_COUNT_CHECK.out.outcome
        paired_reads         // channel: tuple val(meta), path(reads), path(paired_reads): FASTP_TRIMD.out.reads.map
        taxonomy             // channel: tuple val(meta), path(taxonomy): DETERMINE_TAXA_ID.out.taxonomy
        mlst_db              // MLST DB to use with torstens MLST program
        do_srst2_mlst        //
        run_type             // either "update" or "original"

    main:
        ch_versions = Channel.empty() // Used to collect the software versions


        // Creating channel to ensure ID is paired with matching trimmed assembly
        if (run_type=="original") {
            mlst_ch = trimmed_assembly.map{meta, fasta         -> [[id:meta.id], fasta]}\
            .join(scaffold_count_check.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [[id:meta.id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0])\
            .join(taxonomy.map{            meta, taxonomy      -> [[id:meta.id], taxonomy]}, by: [0]).combine(mlst_db)
        } else if (run_type=="update") {
            mlst_ch = trimmed_assembly.map{meta, fasta         -> [[id:meta.id, project_id:meta.project_id], fasta]}\
            .join(scaffold_count_check.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [[id:meta.id, project_id:meta.project_id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])\
            .join(taxonomy.map{             meta, taxonomy      -> [[id:meta.id, project_id:meta.project_id], taxonomy]}, by: [[0][0],[0][1]]).combine(mlst_db)
        }

        // Running standard mlst tool (torstens) on assembly file using provided mlst database location for scemes, profiles, and allele definitions
        MLST (
            mlst_ch
        )
        ch_versions = ch_versions.mix(MLST.out.versions)

        // Creating a channel to pair up the tsv output to the matching taxonomy file, linked on metadata ID
        if (run_type=="original") {
            check_main_mlst_ch = MLST.out.tsv.map{meta, tsv      -> [[id:meta.id], tsv]}\
            .join(taxonomy.map{                   meta, taxonomy -> [[id:meta.id], taxonomy]}, by: [0]).combine(mlst_db)
        } else if (run_type=="update") {
            check_main_mlst_ch = MLST.out.tsv.map{meta, tsv      -> [[id:meta.id, project_id:meta.project_id], tsv]}\
            .join(taxonomy.map{                   meta, taxonomy -> [[id:meta.id, project_id:meta.project_id], taxonomy]}, by: [[0][0],[0][1]]).combine(mlst_db)
        }

        // Checks to see if multiple schemes were found in the sample. Will create _combined.tsv with one ST profile found per line
        CHECK_MLST (
            check_main_mlst_ch
        )
        ch_versions = ch_versions.mix(CHECK_MLST.out.versions)

        // if cdc_phoenix, we want to do srst2 if the standard mlst has reported out as novel_allele, otherwise just skip. 
        // Novel_Profile will not trigger this anymore as no real-time mlst database updates are performed in the CHECK_MLST_WITH_SRST2 process
        if (do_srst2_mlst == true) {
            if (run_type=="original") {
                pre_GET_MLST_SRST2_ch = taxonomy.map                {meta, taxonomy -> [[id:meta.id], taxonomy]}\
                .join(CHECK_MLST.out.status.splitCsv(strip:true).map{meta, status   -> [[id:meta.id], status]},  by: [0]).combine(mlst_db)
            } else if (run_type=="update") {
                pre_GET_MLST_SRST2_ch = taxonomy.map                {meta, taxonomy -> [[id:meta.id, project_id:meta.project_id], taxonomy]}\
                .join(CHECK_MLST.out.status.splitCsv(strip:true).map{meta, status   -> [[id:meta.id, project_id:meta.project_id], status]},  by: [[0][0],[0][1]]).combine(mlst_db)
            }
            // Runs the getMLST portion of the srst2 mlst script to create the required file for srst2 mlst to run. Defines the filesnames and paths to look for in the given database path
            GET_MLST_SRST2 (
                pre_GET_MLST_SRST2_ch
            )
            ch_versions = ch_versions.mix(GET_MLST_SRST2.out.versions)

            // Takes the check statuses from get_MLST step on whether this step has the required files to actually proceed, or should proceed
            if (run_type=="original") {
                mid_srst2_ch = paired_reads.map{                     meta, reads    -> [[id:meta.id], reads]}\
                .join(GET_MLST_SRST2.out.getMLSTs_checker.map{       meta, getMLSTs -> [[id:meta.id], getMLSTs]}, by: [0])\
                .join(GET_MLST_SRST2.out.fastas_checker.map{         meta, fastas   -> [[id:meta.id], fastas]},   by: [0])\
                .join(GET_MLST_SRST2.out.profiles_checker.map{       meta, profiles -> [[id:meta.id], profiles]}, by: [0])\
                .join(CHECK_MLST.out.status.splitCsv(strip:true).map{meta, status   -> [[id:meta.id], status]},   by: [0])
            } else if (run_type=="update") {
                mid_srst2_ch = paired_reads.map{                     meta, reads    -> [[id:meta.id, project_id:meta.project_id], reads]}\
                .join(GET_MLST_SRST2.out.getMLSTs_checker.map{       meta, getMLSTs -> [[id:meta.id, project_id:meta.project_id], getMLSTs]}, by: [[0][0],[0][1]])\
                .join(GET_MLST_SRST2.out.fastas_checker.map{         meta, fastas   -> [[id:meta.id, project_id:meta.project_id], fastas]},   by: [[0][0],[0][1]])\
                .join(GET_MLST_SRST2.out.profiles_checker.map{       meta, profiles -> [[id:meta.id, project_id:meta.project_id], profiles]}, by: [[0][0],[0][1]])\
                .join(CHECK_MLST.out.status.splitCsv(strip:true).map{meta, status   -> [[id:meta.id, project_id:meta.project_id], status]},   by: [[0][0],[0][1]])
            }

            // Identifying mlst genes in trimmed reads
            SRST2_MLST (
                mid_srst2_ch
            )
            ch_versions = ch_versions.mix(SRST2_MLST.out.versions)

            if (run_type=="original") {
                combined_mlst_ch = MLST.out.tsv.map{                        meta, tsv               -> [[id:meta.id], tsv]}\
                .join(SRST2_MLST.out.mlst_results_temp.map{                 meta, mlst_results_temp -> [[id:meta.id], mlst_results_temp]}, by: [0])\
                .join(taxonomy.map{                                         meta, taxonomy          -> [[id:meta.id], taxonomy]},          by: [0])\
                .join(SRST2_MLST.out.empty_checker.splitCsv(strip:true).map{meta, empty_checker     -> [[id:meta.id], empty_checker]},     by: [0]).combine(mlst_db)
            } else if (run_type=="update") {
                combined_mlst_ch = MLST.out.tsv.map{                        meta, tsv               -> [[id:meta.id, project_id:meta.project_id], tsv]}\
                .join(SRST2_MLST.out.mlst_results_temp.map{                 meta, mlst_results_temp -> [[id:meta.id, project_id:meta.project_id], mlst_results_temp]}, by: [[0][0],[0][1]])\
                .join(taxonomy.map{                                         meta, taxonomy          -> [[id:meta.id, project_id:meta.project_id], taxonomy]},          by: [[0][0],[0][1]])\
                .join(SRST2_MLST.out.empty_checker.splitCsv(strip:true).map{meta, empty_checker     -> [[id:meta.id, project_id:meta.project_id], empty_checker]},     by: [[0][0],[0][1]]).combine(mlst_db)
            }

            // Checks to see if multiple schemes were found in the sample. Will create _combined.tsv with one ST profile found per line. Will consolidate if srst2 and standard outputs report the same Types
            CHECK_MLST_WITH_SRST2 (
                combined_mlst_ch
            )
            ch_versions = ch_versions.mix(CHECK_MLST_WITH_SRST2.out.versions)
            if (run_type=="original") {
                checked_mlst_ch = CHECK_MLST_WITH_SRST2.out.checked_MLSTs.map{meta, checked_MLSTs -> [[id:meta.id], checked_MLSTs]}
            } else if (run_type=="update") {
                checked_mlst_ch = CHECK_MLST_WITH_SRST2.out.checked_MLSTs.map{meta, checked_MLSTs -> [[id:meta.id, project_id:meta.project_id], checked_MLSTs]}
            }

        } else {
            if (run_type=="original") {
                checked_mlst_ch = CHECK_MLST.out.checked_MLSTs.map{meta, checked_MLSTs -> [[id:meta.id], checked_MLSTs]}
            } else if (run_type=="update") {
                checked_mlst_ch = CHECK_MLST.out.checked_MLSTs.map{meta, checked_MLSTs -> [[id:meta.id, project_id:meta.project_id], checked_MLSTs]}
            }
        }

    emit:
        checked_MLSTs = checked_mlst_ch
        versions      = ch_versions // channel: [ versions.yml ]
}