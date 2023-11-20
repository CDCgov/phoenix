version 1.0

task combine_phoenix_run {
  input {
    Array[File]? phoenix_tsv_summaries
    Array[File]? griphin_xlsx_summaries
    String? combined_phoenix_tsv_summary_name = "Phoenix_Summary.tsv"
    String? combined_griphin_xlsx_summary_name = "GRiPHin_Summary.xlsx"
  }
  command <<<
    VERSION="v2.1.0-dev"
    echo $VERISON | tee VERSION
    date | tee DATE

    #download phoenix code to get the script from
    nextflow clone cdcgov/phoenix -r $VERSION ./$VERSION/

    #if phoenix tsv files were passed then combine them
    busco_array=()
    if [ ! -z ~{phoenix_tsv_summaries} ]; then
      COUNTER=1
      PHX_ARRAY=(~{sep=',' phoenix_tsv_summaries})
      for i in ${PHX_ARRAY//,/ }; do
        echo "found $i copying to Phoenix_Summary_$COUNTER.tsv"
        cp $i ./Phoenix_Summary_$COUNTER.tsv ;
        #check if this the phoenix summaries were run with CDC_PHOENIX or PHOENIX
        busco_check=$(head -n 1 Phoenix_Summary_$COUNTER.tsv | cut -d$'\t' -f9)
        if [ "$busco_check" == "BUSCO" ]; then
          busco_array+=(true)
          cdc_phoenix="--busco"
        else
          busco_array+=(false)
          cdc_phoenix=""
        fi
        COUNTER=$((COUNTER + 1))
      done
        # Check if all elements have the same boolean value
        # printf prints each element of the array on a new line, then sorts and counts the unique lines using.
        if [[ $(printf "%s\n" "${busco_array[@]}" | sort -u | wc -l) -eq 1 ]]; then
          echo "Values are the same."
          # here the variable cdc_phoenix is the same as the busco argument
          python3 ./$VERSION/bin/Create_phoenix_summary_tsv.py --out ~{combined_phoenix_tsv_summary_name} $cdc_phoenix
        else
          echo "ERROR: Phoenix_Summary.tsv files are a mix of CDC_PHOENIX and PHOENIX outputs and they need to be the same."
          exit 1
        fi

      #check if the file is empty (aka has something in the 2nd line) and if it is then delete it to cause failure
      if [ "$(wc -l <~{combined_phoenix_tsv_summary_name})" -eq 1 ]; then
        echo "file only contains a single line"
        rm -r ~{combined_phoenix_tsv_summary_name}
        exit 1
      fi
    # if array is empty
    else
      echo "WARNING: No Phoenix_Summary.tsv files provided skipping Phoenix_Summary.tsv combining step."
    fi

    #if griphin xlsx files were passed then combine them
    if [ ! -z ~{griphin_xlsx_summaries} ]; then
      COUNTER=1
      GRIPHIN_ARRAY=(~{sep=',' griphin_xlsx_summaries})
      for i in ${GRIPHIN_ARRAY//,/ }; do
        echo "found $i copying to GRiPHin_Summary_$COUNTER.xlsx"
        cp $i ./GRiPHin_Summary_$COUNTER.xlsx ;
        COUNTER=$((COUNTER + 1))
      done

      ## combine griphin reports. in the script it determines if phx or cdc_phx was run.
      python3 ./$VERSION/bin/terra_combine_griphin.py --out ~{combined_griphin_xlsx_summary_name}

    # if array is empty
    else
      echo "WARNING: No Phoenix_Summary.tsv files provided skipping Phoenix_Summary.tsv combining step."
    fi

  >>>
  output {
    File?   phoenix_tsv_summary  = "~{combined_phoenix_tsv_summary_name}"
    File?   griphin_xlsx_summary = "~{combined_griphin_xlsx_summary_name}"
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
