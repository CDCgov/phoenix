version 1.0

import "../tasks/combine_phoenix_run.wdl" as combine_phoenix_run_nf

workflow combine_phoenix_output {
  meta {
    description: "A WDL wrapper to combine the output of a phoenix run."
  }
  input {
    Array[File] phoenix_tsv_summary
    String? output_file
    String? cdc
  }
  call combine_phoenix_run_nf.combine_phoenix_run {
    input:
      phoenix_tsv_summary = phoenix_tsv_summary,
      output_file         = output_file,
      cdc                 = cdc
  }
  output {
    #phoenix summary output
    File   phoenix_tsv_summary = combine_phoenix_run.phoenix_tsv_summary
    String phoenix_version     = combine_phoenix_run.phoenix_version
    String phoenix_docker      = combine_phoenix_run.phoenix_docker
    String analysis_date       = combine_phoenix_run.analysis_date
  }
}