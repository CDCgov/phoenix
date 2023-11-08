version 1.0

import "../tasks/combine_phoenix_run.wdl" as combine_phoenix_run

workflow combine_phoenix_output {
  meta {
    description: "A WDL wrapper to combine the output of a phoenix run."
  }
  input {
    Array[File] phoenix_tsv_summary
    String? output_file
    String? cdc
  }
  call combine_phoenix_run.combine_phoenix {
    input:
      phoenix_tsv_summary = phoenix_tsv_summary,
      output_file         = output_file,
      cdc                 = cdc
  }
  output {
    #phoenix summary output
    File   phoenix_tsv_summary = combine_phoenix.phoenix_tsv_summary
    String phoenix_version     = combine_phoenix.phoenix_version
    String phoenix_docker      = combine_phoenix.phoenix_docker
    String analysis_date       = combine_phoenix.analysis_date
  }
}