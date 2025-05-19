version 1.0

import "../tasks/update_phoenix.wdl" as update_phoenix_nf

workflow update_phoenix {
  meta {
    description: "A WDL wrapper to update the AR gene calls and MLST information for isolates."
  }
  input {
    Array[File]? phoenix_tsv_summaries
    Array[File]? griphin_xlsx_summaries
    Array[File]? griphin_tsv_summaries
    Array[File]? ncbi_biosample_excel_files
    Array[File]? ncbi_sra_excel_files
    String? combined_phoenix_tsv_prefix
    String? combined_griphin_xlsx_prefix
    String? combined_griphin_tsv_prefix
    String? combined_ncbi_biosample_xlsx_prefix
    String? combined_ncbi_sra_xlsx_prefix
  }
  call update_phoenix_nf.update_phoenix {
    input:
      phoenix_tsv_summaries               = phoenix_tsv_summaries,
      griphin_xlsx_summaries              = griphin_xlsx_summaries,
      griphin_tsv_summaries               = griphin_tsv_summaries,
      ncbi_biosample_excel_files          = ncbi_biosample_excel_files,
      ncbi_sra_excel_files                = ncbi_sra_excel_files,
      combined_phoenix_tsv_prefix         = combined_phoenix_tsv_prefix,
      combined_griphin_xlsx_prefix        = combined_griphin_xlsx_prefix,
      combined_griphin_tsv_prefix         = combined_griphin_tsv_prefix,
      combined_ncbi_biosample_xlsx_prefix = combined_ncbi_biosample_xlsx_prefix,
      combined_ncbi_sra_xlsx_prefix       = combined_ncbi_sra_xlsx_prefix
  }
  output {
    #phoenix summary output
    File?  phoenix_tsv_summary     = update_phoenix.phoenix_tsv_summary
    File?  griphin_xlsx_summary    = update_phoenix.griphin_xlsx_summary
    File?  griphin_tsv_summary     = update_phoenix.griphin_tsv_summary
    File?  biosample_excel_summary = update_phoenix.biosample_excel_summary
    File?  sra_excel_summary       = update_phoenix.sra_excel_summary
    String phoenix_version         = update_phoenix.phoenix_version
    String phoenix_docker          = update_phoenix.phoenix_docker
    String analysis_date           = update_phoenix.analysis_date
  }
}