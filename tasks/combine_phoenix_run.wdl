version 1.0

task combine_phoenix_run {
  input {
    Array[File]? phoenix_tsv_summaries
    Array[File]? griphin_summaries
    String? phoenix_tsv_summary_name = "Phoenix_Summary.tsv"
    String? griphin_xlsx_name = "GRiPHin_Summary.xlsx"
  }
  command <<<
    VERSION="v2.1.0-dev"
    echo $VERISON | tee VERSION
    date | tee DATE

    #download phoenix code to get the script from
    nextflow clone cdcgov/phoenix -r $VERSION ./$VERSION/

    #if phoenix tsv files were passed then combine them
    if [ ! -z ~{phoenix_tsv_summaries} ]; then
      COUNTER=1
      ARRAY=(~{sep=',' phoenix_tsv_summaries})
      for i in ${ARRAY//,/ }; do
        echo "found $i copying to Phoenix_Summary_$COUNTER.tsv"
        cp $i ./Phoenix_Summary_$COUNTER.tsv ;
        COUNTER=$((COUNTER + 1))
        #check if this the phoenix summaries were run with CDC_PHOENIX or PHOENIX

      done

      ## here ~{cdc} is the same as the busco argument
      python3 ./$VERSION/bin/Create_phoenix_summary_tsv.py --out ~{phoenix_summary_name} $cdc_phoenix

      #check if the file is empty (aka has something in the 2nd line) and if it is then delete it to cause failure
      if [ "$(wc -l <~{phoenix_summary_name})" -eq 1 ]; then
        echo "file only contains a single line"
        rm -r ~{phoenix_summary_name}
        exit 1
      fi
    fi

    #if griphin xlsx files were passed then combine them
    if [ ! -z ~{griphin_summaries} ]; then
      COUNTER=1
      ARRAY=(~{sep=',' griphin_summaries})
      for i in ${ARRAY//,/ }; do
        echo "found $i copying to GRiPHin_Summary_$COUNTER.xlsx"
        cp $i ./GRiPHin_Summary_$COUNTER.xlsx ;
        COUNTER=$((COUNTER + 1))
      done

      ## here ~{cdc} is the same as the busco argument
      python3 ./$VERSION/bin/terra_combine_griphin.py --out ~{griphin_xlsx_name} ~{cdc}
    fi
  
  >>>
  output {
    File?   phoenix_tsv_summary  = "~{phoenix_summary_name}"
    File?   griphin_xlsx_summary = "~{griphin_xlsx_name}"
    String  phoenix_version      = read_string("VERSION")
    String  phoenix_docker       = "quay.io/jvhagey/phoenix:2.0.2"
    String  analysis_date        = read_string("DATE")
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
