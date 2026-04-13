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
        trimmed_assembly     // channel: tuple val(meta), path(assembly)
        scaffold_count_check // SCAFFOLD_COUNT_CHECK.out.outcome
        paired_reads         // channel: tuple val(meta), path(reads)
        taxonomy             // channel: tuple val(meta), path(taxonomy)
        mlst_db              // MLST DB
        do_srst2_mlst        // Boolean
        run_type             // "update" or "original"

    main:
        ch_versions = Channel.empty()

        // 1. Prepare input for standard MLST
        if (run_type == "original") {
            mlst_ch = trimmed_assembly.map{ meta, fasta -> [[id:meta.id], fasta] }
                .join(scaffold_count_check.splitCsv(strip:true, by:5).map{ meta, f -> [[id:meta.id], [f[0][0], f[1][0], f[2][0], f[3][0], f[4][0]]] }, by: [0])
                .filter { meta, fasta, f -> f[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering." }
                .join(taxonomy.map{ meta, tax -> [[id:meta.id], tax] }, by: [0])
                .map{ meta, fasta, f, tax -> [meta, fasta, tax] }
        } else {
            mlst_ch = trimmed_assembly.map{ meta, fasta -> [[id:meta.id, project_id:meta.project_id], fasta] }
                .join(scaffold_count_check.map{ meta, outcome ->
                    def content = file(outcome.toString()).text
                    def passed = !content.contains("FAILED")
                    [[id:meta.id, project_id:meta.project_id], [outcome, passed]]
                }, by: [[0][0],[0][1]])
                .filter{ meta, fasta, data -> data[1] == true }
                .map{ meta, fasta, data -> [meta, fasta] }
                .join(taxonomy.map{ meta, tax -> [[id:meta.id, project_id:meta.project_id], tax] }, by: [[0][0],[0][1]])
        }

        // 2. Run Standard MLST
        MLST ( mlst_ch )
        ch_versions = ch_versions.mix(MLST.out.versions)

        // 3. Run Standard Checker
        if (run_type == "original") {
            check_main_mlst_ch = MLST.out.tsv.map{ meta, tsv -> [[id:meta.id], tsv] }
                .join(taxonomy.map{ meta, tax -> [[id:meta.id], tax] }, by: [0]).combine(mlst_db)
        } else {
            check_main_mlst_ch = MLST.out.tsv.map{ meta, tsv -> [[id:meta.id, project_id:meta.project_id], tsv] }
                .join(taxonomy.map{ meta, tax -> [[id:meta.id, project_id:meta.project_id], tax] }, by: [[0][0],[0][1]]).combine(mlst_db)
        }

        CHECK_MLST ( check_main_mlst_ch )
        ch_versions = ch_versions.mix(CHECK_MLST.out.versions)

        // 4. Handle SRST2 Logic with Fallback
        if (do_srst2_mlst == true) {
            if (run_type == "original") {
                pre_get_ch = taxonomy.map{ m, t -> [[id:m.id], t] }.join(CHECK_MLST.out.status.splitCsv(strip:true).map{ m, s -> [[id:m.id], s] }, by: [0]).combine(mlst_db)
            } else {
                pre_get_ch = taxonomy.map{ m, t -> [[id:m.id, project_id:m.project_id], t] }.join(CHECK_MLST.out.status.splitCsv(strip:true).map{ m, s -> [[id:m.id, project_id:m.project_id], s] }, by: [[0][0],[0][1]]).combine(mlst_db)
            }

            GET_MLST_SRST2 ( pre_get_ch )
            
            if (run_type == "original") {
                mid_srst2_ch = paired_reads.map{ m, r -> [[id:m.id], r] }
                    .join(GET_MLST_SRST2.out.getMLSTs_checker.map{ m, f -> [[id:m.id], f] }, by: [0])
                    .join(GET_MLST_SRST2.out.fastas_checker.map{ m, f -> [[id:m.id], f] }, by: [0])
                    .join(GET_MLST_SRST2.out.profiles_checker.map{ m, f -> [[id:m.id], f] }, by: [0])
                    .join(CHECK_MLST.out.status.splitCsv(strip:true).map{ m, s -> [[id:m.id], s] }, by: [0])
            } else {
                mid_srst2_ch = paired_reads.map{ m, r -> [[id:m.id, project_id:m.project_id], r] }
                    .join(GET_MLST_SRST2.out.getMLSTs_checker.map{ m, f -> [[id:m.id, project_id:m.project_id], f] }, by: [[0][0],[0][1]])
                    .join(GET_MLST_SRST2.out.fastas_checker.map{ m, f -> [[id:m.id, project_id:m.project_id], f] }, by: [[0][0],[0][1]])
                    .join(GET_MLST_SRST2.out.profiles_checker.map{ m, f -> [[id:m.id, project_id:m.project_id], f] }, by: [[0][0],[0][1]])
                    .join(CHECK_MLST.out.status.splitCsv(strip:true).map{ m, s -> [[id:m.id, project_id:m.project_id], s] }, by: [[0][0],[0][1]])
            }

            SRST2_MLST ( mid_srst2_ch )

            if (run_type == "original") {
                combined_mlst_ch = MLST.out.tsv.map{ m, t -> [[id:m.id], t] }
                    .join(SRST2_MLST.out.mlst_results_temp.map{ m, f -> [[id:m.id], f] }, by: [0])
                    .join(taxonomy.map{ m, t -> [[id:m.id], t] }, by: [0])
                    .join(SRST2_MLST.out.empty_checker.splitCsv(strip:true).map{ m, e -> [[id:m.id], e] }, by: [0]).combine(mlst_db)
            } else {
                combined_mlst_ch = MLST.out.tsv.map{ m, t -> [[id:m.id, project_id:m.project_id], t] }
                    .join(SRST2_MLST.out.mlst_results_temp.map{ m, f -> [[id:m.id, project_id:m.project_id], f] }, by: [[0][0],[0][1]])
                    .join(taxonomy.map{ m, t -> [[id:m.id, project_id:m.project_id], t] }, by: [[0][0],[0][1]])
                    .join(SRST2_MLST.out.empty_checker.splitCsv(strip:true).map{ m, e -> [[id:m.id, project_id:m.project_id], e] }, by: [[0][0],[0][1]]).combine(mlst_db)
            }

            CHECK_MLST_WITH_SRST2 ( combined_mlst_ch )
            ch_versions = ch_versions.mix(CHECK_MLST_WITH_SRST2.out.versions)

            // --- THE BULLETPROOF FALLBACK ---
            def assembly_only = CHECK_MLST.out.checked_MLSTs.map{ meta, f -> [meta.id, meta, f] }
            def srst2_checked = CHECK_MLST_WITH_SRST2.out.checked_MLSTs.map{ meta, f -> [meta.id, meta, f] }

            checked_mlst_ch = assembly_only.join(srst2_checked, remainder: true)
                .map { it ->
                    // it[0] is the ID
                    // it[1] is meta_asm, it[2] is file_asm
                    // it[3] is meta_srst2, it[4] is file_srst2
                    
                    def final_meta = it.size() > 3 && it[3] != null ? it[3] : it[1]
                    def final_file = it.size() > 4 && it[4] != null ? it[4] : it[2]
                    
                    return [ final_meta, final_file ]
                }

        } else {
            // Standard Path if do_srst2 is false
            checked_mlst_ch = CHECK_MLST.out.checked_MLSTs
        }

    emit:
        checked_MLSTs = checked_mlst_ch
        //checked_MLSTs = checked_mlst_ch.view { "DEBUG DO_MLST EMIT: $it" }
        versions      = ch_versions
}