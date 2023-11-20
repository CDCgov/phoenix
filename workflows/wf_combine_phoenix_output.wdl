version 1.0

import "../tasks/combine_phoenix_run.wdl" as combine_phoenix_run_nf

workflow combine_phoenix_output {
  meta {
    description: "A WDL wrapper to combine the output of a phoenix run."
  }
  input {
    Array[File]? phoenix_tsv_summaries
    Array[File]? griphin_xlsx_summaries
    String? phoenix_tsv_summary_name
    String? griphin_xlsx_name
  }
  call combine_phoenix_run_nf.combine_phoenix_run {
    input:
      phoenix_tsv_summaries    = phoenix_tsv_summaries,
      griphin_xlsx_summaries   = griphin_xlsx_summaries,
      phoenix_tsv_summary_name = phoenix_tsv_summary_name,
      griphin_xlxs_name        = griphin_xlsx_name
  }
  output {
    #phoenix summary output
    File?  phoenix_tsv_summary  = combine_phoenix_run.phoenix_tsv_summary
    File?  griphin_xlsx_summary = combine_phoenix_run.griphin_xlsx_summary
    String phoenix_version      = combine_phoenix_run.phoenix_version
    String phoenix_docker       = combine_phoenix_run.phoenix_docker
    String analysis_date        = combine_phoenix_run.analysis_date
  }
}