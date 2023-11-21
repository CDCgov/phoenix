version 1.0

task combine_phoenix_run {
  input {
    Array[File]? phoenix_tsv_summaries
<<<<<<< HEAD
    Array[File]? griphin_xlsx_summaries
    Array[File]? griphin_tsv_summaries
    String? combined_phoenix_tsv_summary_name = "Phoenix_Summary.tsv"
    String? combined_griphin_xlsx_summary_name = "GRiPHin_Summary.xlsx"
    String? combined_griphin_tsv_summary_name = "GRiPHin_Summary.tsv"
=======
    Array[File]? griphin_summaries
    String? phoenix_tsv_summary_name = "Phoenix_Summary.tsv"
    String? griphin_xlsx_name = "GRiPHin_Summary.xlsx"
>>>>>>> 72945a7dfeb7e68e21df20aa1f03218b6467181d
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
      python3 ./$VERSION/bin/Create_phoenix_summary_tsv.py --out ~{phoenix_tsv_summary_name} $cdc_phoenix

      #check if the file is empty (aka has something in the 2nd line) and if it is then delete it to cause failure
      if [ "$(wc -l <~{phoenix_tsv_summary_name})" -eq 1 ]; then
        echo "file only contains a single line"
        rm -r ~{phoenix_tsv_summary_name}
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

      ## combine griphin summaries. In the script it determines if phx or cdc_phx was run.
      python3 ./$VERSION/bin/terra_combine_griphin.py --out ~{combined_griphin_xlsx_summary_name}

      # If GRiPHin files were passed, but not a summary made at the end then throw an error
      if [ ! -s "~{combined_griphin_xlsx_summary_name}" ] && [ ! -f "~{combined_griphin_xlsx_summary_name}" ]; then
        echo "ERROR: GRiPHin excel files were passed, but no combination file was made."
        exit 1
      fi
    # if array is empty
    else
      echo "WARNING: No GRiPHin_Summary.xlsx files provided skipping GRiPHin_Summary.xlsx combining step."
    fi

    #if griphin tsv files were passed then combine them
    busco_gripin_array=()
    if [ ! -z "~{sep=',' griphin_tsv_summaries}" ]; then
      COUNTER=1
      GRIPHIN_ARRAY_TSV=(~{sep=',' griphin_tsv_summaries})
      for i in ${GRIPHIN_ARRAY_TSV//,/ }; do
        echo "found $i copying to GRiPHin_${COUNTER}_Summary.tsv"
        cp $i ./GRiPHin_${COUNTER}_Summary.tsv ;
        COUNTER=$((COUNTER + 1))
      done

      ## combine griphin reports. In the script it determines if phx or cdc_phx was run.
      python3 ./$VERSION/bin/terra_combine_griphin_tsv.py --out ~{combined_griphin_tsv_summary_name}

      # If GRiPHin files were passed, but not a summary made at the end then throw an error
      if [ ! -s "~{combined_griphin_tsv_summary_name}" ] && [ ! -f "~{combined_griphin_tsv_summary_name}" ]; then
        echo "ERROR: GRiPHin tsv files were passed, but no combination file was made."
        exit 1
      fi
    # if array is empty
    else
      echo "WARNING: No GRiPHin_Summary.tsv files provided skipping GRiPHin_Summary.tsv combining step."
    fi

  >>>
  output {
    File?   phoenix_tsv_summary  = "~{combined_phoenix_tsv_summary_name}"
    File?   griphin_xlsx_summary = "~{combined_griphin_xlsx_summary_name}"
    File?   griphin_tsv_summary  = "~{combined_griphin_tsv_summary_name}"
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
