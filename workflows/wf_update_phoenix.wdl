version 1.0

import "../tasks/update_phoenix.wdl" as update_phoenix_nf

workflow update_phoenix {
  meta {
    description: "A WDL wrapper to update the AR gene calls and MLST information for isolates."
  }
  input {
    String  samplename
    String  project_directory
    Int?    coverage
  }
  call update_phoenix_nf.update_phoenix {
    input:
      samplename        = samplename,
      coverage          = coverage,
      project_directory = project_directory
  }
  output {
    #phoenix summary output values
    File?   work_files                        = update_phoenix.work_files
    String  phoenix_version                   = update_phoenix.phoenix_version
    String  phoenix_docker                    = update_phoenix.phoenix_docker
    String  analysis_date                     = update_phoenix.analysis_date
    String  qc_outcome                        = update_phoenix.qc_outcome
    String  warnings                          = update_phoenix.warnings
    String  final_taxa_id                     = update_phoenix.final_taxa_id
    String  taxa_source                       = update_phoenix.taxa_source
    String? shigapass_taxa                    = update_phoenix.shigapass_taxa
    String  mlst_scheme_1                     = update_phoenix.mlst_scheme_1
    String  mlst_1                            = update_phoenix.mlst_1
    String  mlst1_ncbi                        = update_phoenix.mlst1_ncbi
    String  mlst_scheme_2                     = update_phoenix.mlst_scheme_2
    String  mlst_2                            = update_phoenix.mlst_2
    String  mlst2_ncbi                        = update_phoenix.mlst2_ncbi
    String  gamma_beta_lactam_genes           = update_phoenix.gamma_beta_lactam_genes
    String  other_ar_genes                    = update_phoenix.other_ar_genes
    String  amrfinder_point_mutations         = update_phoenix.amrfinder_point_mutations
    String  amrfinder_amr_classes             = update_phoenix.amrfinder_amr_classes
    String  amrfinder_amr_subclasses          = update_phoenix.amrfinder_amr_subclasses
    String  amrfinder_core_genes              = update_phoenix.amrfinder_core_genes
    String  amrfinder_plus_genes              = update_phoenix.amrfinder_plus_genes
    String  amrfinder_stress_genes            = update_phoenix.amrfinder_stress_genes
    String  amrfinder_virulence_genes         = update_phoenix.amrfinder_virulence_genes
    String  amrfinder_beta_lactam_genes       = update_phoenix.amrfinder_beta_lactam_genes
    String  hypervirulence_genes              = update_phoenix.hypervirulence_genes
    String  plasmid_incompatibility_replicons = update_phoenix.plasmid_incompatibility_replicons
    String  qc_issues                         = update_phoenix.qc_issues
    #phoenix ani
    File?  reformated_fast_ani      = update_phoenix.reformated_fast_ani
    #phoenix quast and mlst
    File   mlst_tsv                 = update_phoenix.mlst_tsv
    # cdc_phoenix busco and srst2 - optional for PHOENIX, SCAFFOLDS and SRA entries
    File?  srst2                    = update_phoenix.srst2
    #phoenix gamma
    File   gamma_ar_calls           = update_phoenix.gamma_ar_calls
    File   blat_ar_calls            = update_phoenix.blat_ar_calls
    #phoenix output
    File   summary_line             = update_phoenix.summary_line
    File   synopsis                 = update_phoenix.synopsis
    File?  best_taxa_id             = update_phoenix.best_taxa_id
    #phoenix amrfinder
    File   amrfinder_mutations      = update_phoenix.amrfinder_mutations
    File?  amrfinder_taxa_match     = update_phoenix.amrfinder_taxa_match
    File?  amrfinder_hits           = update_phoenix.amrfinder_hits
    #species specific
    File?  shigapass_summary       = update_phoenix.shigapass_summary
    #phoenix summary output
    File   updater_log             = update_phoenix.updater_log
    File   phoenix_tsv_summary     = update_phoenix.phoenix_tsv_summary
    File   griphin_excel_summary   = update_phoenix.griphin_excel_summary
    File   griphin_tsv_summary     = update_phoenix.griphin_tsv_summary
    String phoenix_version         = update_phoenix.phoenix_version
    String phoenix_docker          = update_phoenix.phoenix_docker
    String analysis_date           = update_phoenix.analysis_date
  }
}