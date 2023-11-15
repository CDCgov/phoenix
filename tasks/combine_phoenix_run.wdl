version 1.0

task combine_phoenix_run {
  input {
    Array[File] phoenix_tsv_summaries
    String? output_file = "Phoenix_Summary.tsv"
    String? cdc = false
  }
  command <<<
    VERSION="v2.1.0-dev"
    echo $VERISON | tee VERSION
    date | tee DATE

    #download phoenix code to get the script from
    nextflow clone cdcgov/phoenix -r $VERSION ./$VERSION/

    COUNTER=1
    ARRAY=(~{sep=',' phoenix_tsv_summaries})
    for i in ${ARRAY//,/ }; do
      echo "found $i copying to Phoenix_Summary_$COUNTER.tsv"
      cp $i ./Phoenix_Summary_$COUNTER.tsv ;
      COUNTER=$((COUNTER + 1))
    done

    ## here ~{cdc} is the same as the busco argument
    python3 ./$VERSION/bin/Create_phoenix_summary_tsv.py --out ~{output_file} ~{cdc}

    #check if the file is empty (aka has something in the 2nd line) and if it is then delete it to cause failure
    if [ "$(wc -l <~{output_file})" -eq 1 ]; then
      echo "file only contains a single line"
      rm -r ~{output_file}
    fi
  
  >>>
  output {
    File    phoenix_tsv_summary = "~{output_file}"
    String  phoenix_version     = read_string("VERSION")
    String  phoenix_docker      = "quay.io/jvhagey/phoenix:2.0.2"
    String  analysis_date       = read_string("DATE")
  }
  runtime {
    docker: "quay.io/jvhagey/phoenix:2.0.2"
    memory: "8 GB"
    cpu: 1
    disks:  "local-disk 100 SSD"
    maxRetries: 0
    preemptible: 0
  }
}
