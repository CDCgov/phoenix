version 1.0

import "../tasks/combine_phoenix_run.wdl" as combine_phoenix_run_nf

workflow combine_phoenix_output {
  meta {
    description: "A WDL wrapper to combine the output of a phoenix run."
  }
  input {
    Array[File]? phoenix_tsv_summaries
    Array[File]? griphin_xlsx_summaries
    Array[File]? griphin_tsv_summaries
    Array[File]? ncbi_biosample_attributes_excel_files
    Array[File]? ncbi_sra_excel_files
    String? combined_phoenix_tsv_prefix
    String? combined_griphin_xlsx_prefix
    String? combined_griphin_tsv_prefix
    String? combined_ncbi_biosample_xlsx_prefix
    String? combined_sra_biosample_xlsx_prefix
  }
  call combine_phoenix_run_nf.combine_phoenix_run {
    input:
      phoenix_tsv_summaries              = phoenix_tsv_summaries,
      griphin_xlsx_summaries             = griphin_xlsx_summaries,
      griphin_tsv_summaries              = griphin_tsv_summaries,
      ncbi_biosample_attributes_excel_files = ncbi_biosample_attributes_excel_files,
      ncbi_sra_excel_files = ncbi_sra_excel_files,
      combined_phoenix_tsv_prefix  = combined_phoenix_tsv_prefix,
      combined_griphin_xlsx_prefix = combined_griphin_xlsx_prefix,
      combined_griphin_tsv_prefix  = combined_griphin_tsv_prefix,
      combined_ncbi_biosample_xlsx_prefix = combined_ncbi_biosample_xlsx_prefix,
      combined_sra_biosample_xlsx_prefix = combined_sra_biosample_xlsx_prefix
  }
  output {
    #phoenix summary output
    File?  phoenix_tsv_summary     = combine_phoenix_run.phoenix_tsv_summary
    File?  griphin_xlsx_summary    = combine_phoenix_run.griphin_xlsx_summary
    File?  griphin_tsv_summary     = combine_phoenix_run.griphin_tsv_summary
    File?  biosample_excel_summary = combine_phoenix_run.biosample_excel_summary
    File?  sra_excel_summary       = combine_phoenix_run.sra_excel_summary
    String phoenix_version         = combine_phoenix_run.phoenix_version
    String phoenix_docker          = combine_phoenix_run.phoenix_docker
    String analysis_date           = combine_phoenix_run.analysis_date
  }
}