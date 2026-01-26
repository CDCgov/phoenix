version 1.0

task combine_phoenix_run {
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
  command <<<
    version="v2.2.1"
    echo $version | tee VERSION
    date | tee DATE

    #download phoenix code to get the script from
    nextflow clone cdcgov/phoenix -r $version ./$version/

    # create file name depending on the user picked prefix
    if [ ! -z "~{combined_phoenix_tsv_prefix}" ]; then
      out_command="--out ~{combined_phoenix_tsv_prefix}_Phoenix_Summary.tsv"
      combined_phoenix_tsv_summary_name="~{combined_phoenix_tsv_prefix}_Phoenix_Summary.tsv"
    else
      out_command="--out Phoenix_Summary.tsv"
      combined_phoenix_tsv_summary_name="Phoenix_Summary.tsv"
    fi
    if [ ! -z "~{combined_griphin_xlsx_prefix}" ]; then
      out_griphin_xlsx_command="--out ~{combined_griphin_xlsx_prefix}"
      combined_griphin_xlsx_summary_name="~{combined_griphin_xlsx_prefix}_GRiPHin_Summary.xlsx"
    else
      out_griphin_xlsx_command=""
      combined_griphin_xlsx_summary_name="GRiPHin_Summary.xlsx"
    fi
    if [ ! -z "~{combined_griphin_tsv_prefix}" ]; then
      out_griphin_tsv_command="--out ~{combined_griphin_tsv_prefix}"
      combined_griphin_tsv_summary_name="~{combined_griphin_tsv_prefix}_GRiPHin_Summary.tsv"
    else
      out_griphin_tsv_command=""
      combined_griphin_tsv_summary_name="GRiPHin_Summary.tsv"
    fi
    if [ ! -z "~{combined_ncbi_sra_xlsx_prefix}" ]; then
      out_ncbi_sra_command="--sra_output ~{combined_ncbi_sra_xlsx_prefix}"
      combined_ncbi_sra_summary_name="~{combined_ncbi_sra_xlsx_prefix}_Sra_Microbe.1.0.xlsx"
    else
      out_ncbi_sra_command=""
      combined_ncbi_sra_summary_name="Sra_Microbe.1.0.xlsx"
    fi
    if [ ! -z "~{combined_ncbi_biosample_xlsx_prefix}" ]; then
      out_ncbi_biosample_command="--biosample_output ~{combined_ncbi_biosample_xlsx_prefix}"
      combined_ncbi_biosample_summary_name="~{combined_ncbi_biosample_xlsx_prefix}_BiosampleAttributes_Microbe.1.0.xlsx"
    else
      out_ncbi_biosample_command=""
      combined_ncbi_biosample_summary_name="BiosampleAttributes_Microbe.1.0.xlsx"
    fi 

    #if phoenix tsv files were passed then combine them
    busco_array=()
    if [ ! -z "~{sep=',' phoenix_tsv_summaries}" ]; then
      echo "Combining and creating ${combined_phoenix_tsv_summary_name}"
      COUNTER=1
      PHX_ARRAY=(~{sep=',' phoenix_tsv_summaries})
      for i in ${PHX_ARRAY//,/ }; do
        echo "found $i copying to Phoenix_Summary_${COUNTER}.tsv"
        cp $i ./Phoenix_Summary_${COUNTER}.tsv ;
        #check if the phoenix summaries were run with CDC_PHOENIX or PHOENIX. They need to be the same.
        busco_check=$(head -n 1 Phoenix_Summary_${COUNTER}.tsv | cut -d$'\t' -f9)
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
          echo "Phoenix_Summary.tsv files passed check for the same entry point. Starting to combine files."
          # here the variable cdc_phoenix is the same as the busco argument
          python3 ./$version/bin/Create_phoenix_summary_tsv.py $out_command $cdc_phoenix
        else
          echo "ERROR: Phoenix_Summary.tsv files are a mix of CDC_PHOENIX and PHOENIX outputs and they need to be the same."
          exit 1
        fi

      #check if the file is empty (aka has something in the 2nd line) and if it is then delete it to cause failure
      if [[ "$(wc -l <${combined_phoenix_tsv_summary_name})" -eq 1 ]]; then
        echo "ERROR: Phoenix_Summary.tsv only contains a single line. Combination failed."
        rm -r ${combined_phoenix_tsv_summary_name}
        exit 1
      fi
    # if array is empty
    else
      echo "WARNING: No Phoenix_Summary.tsv files provided skipping Phoenix_Summary.tsv combining step."
    fi

    #if griphin xlsx files were passed then combine them
    if [ ! -z "~{sep=',' griphin_xlsx_summaries}" ]; then
      echo "Combining and creating ${combined_griphin_xlsx_summary_name}"
      COUNTER=1
      GRIPHIN_ARRAY=(~{sep=',' griphin_xlsx_summaries})
      for i in ${GRIPHIN_ARRAY//,/ }; do
        echo "found $i copying to GRiPHin_${COUNTER}_Summary.xlsx"
        cp $i ./GRiPHin_${COUNTER}_Summary.xlsx ;
        COUNTER=$((COUNTER + 1))
      done

      ## combine griphin summaries. In the script it determines if phx or cdc_phx was run.
      python3 ./$version/bin/terra_combine_griphin.py $out_griphin_xlsx_command

      # If GRiPHin files were passed, but not a summary made at the end then throw an error
      if [ ! -s "${combined_griphin_xlsx_summary_name}" ]; then
        echo "ERROR: GRiPHin excel files were passed, but no combination file was made."
        ls
        exit 1
      fi
    # if array is empty
    else
      echo "WARNING: No GRiPHin_Summary.xlsx files provided skipping GRiPHin_Summary.xlsx combining step."
    fi

    #if griphin tsv files were passed then combine them
    busco_gripin_array=()
    if [ ! -z "~{sep=',' griphin_tsv_summaries}" ]; then
      echo "Combining and creating ${combined_griphin_tsv_summary_name}"
      COUNTER=1
      GRIPHIN_ARRAY_TSV=(~{sep=',' griphin_tsv_summaries})
      for i in ${GRIPHIN_ARRAY_TSV//,/ }; do
        echo "found $i copying to GRiPHin_${COUNTER}_Summary.tsv"
        cp $i ./GRiPHin_${COUNTER}_Summary.tsv ;
        COUNTER=$((COUNTER + 1))
      done

      ## combine griphin reports. In the script it determines if phx or cdc_phx was run.
      python3 ./$version/bin/terra_combine_griphin_tsv.py $out_griphin_tsv_command

      # If GRiPHin files were passed, but not a summary made at the end then throw an error
      if [ ! -s "${combined_griphin_tsv_summary_name}" ]; then
        echo "ERROR: GRiPHin tsv files were passed, but no combination file was made."
        ls
        exit 1
      fi
    # if array is empty
    else
      echo "WARNING: No GRiPHin_Summary.tsv files provided skipping GRiPHin_Summary.tsv combining step."
    fi

  #if ncbi biosampe excel files were passed then combine them
  if [ ! -z "~{sep=',' ncbi_biosample_excel_files}" ]; then
    echo "Combining and creating ${combined_ncbi_biosample_summary_name}"
    COUNTER=1
    BIOSAMPLE_ARRAY_EXCEL=(~{sep=',' ncbi_biosample_excel_files})
    for i in ${BIOSAMPLE_ARRAY_EXCEL//,/ }; do
      echo "found $i copying to BiosampleAttributes_${COUNTER}_Microbe.1.0.xlsx"
      cp $i ./BiosampleAttributes_${COUNTER}_Microbe.1.0.xlsx ;
      COUNTER=$((COUNTER + 1))
    done

    ## combine ncbi biosample excel files. In the script it determines if phx or cdc_phx was run.
    python3 ./$version/bin/terra_combine_ncbi_excel.py $out_ncbi_biosample_command

    # If NCBI biosample files were passed, but not a summary made at the end then throw an error
    if [ ! -s "${combined_ncbi_biosample_summary_name}" ]; then
      echo "ERROR: BiosampleAttributes_Microbe.1.0.xlsx files were passed, but no combination file was made."
      ls
      exit 1
    fi
  # if array is empty
  else
    echo "WARNING: No NCBI excel files provided skipping NCBI excel combining steps."
  fi

  if [ ! -z "~{sep=',' ncbi_sra_excel_files}" ]; then
    echo "Combining and creating ${combined_ncbi_sra_summary_name}"
    COUNTER=1
    SRA_ARRAY_EXCEL=(~{sep=',' ncbi_sra_excel_files})
    for i in ${SRA_ARRAY_EXCEL//,/ }; do
      echo "found $i copying to Sra_${COUNTER}_Microbe.1.0.xlsx"
      cp $i ./Sra_${COUNTER}_Microbe.1.0.xlsx ;
      COUNTER=$((COUNTER + 1))
    done

    ## combine sra excel files. In the script it determines if phx or cdc_phx was run.
    python3 ./$version/bin/terra_combine_ncbi_excel.py $out_ncbi_sra_command

    # If NCBI biosample files were passed, but not a summary made at the end then throw an error
    if [ ! -s "${combined_ncbi_sra_summary_name}" ]; then
      echo "ERROR: Sra_Microbe.1.0.xlsx files were passed, but no combination file was made."
      ls
      exit 1
    fi
  # if array is empty
  else
    echo "WARNING: No NCBI excel files provided skipping NCBI excel combining steps."
  fi


  # series of checks to finish up
  #check at least one file type was passed, if not then fail.
  if [ -z "~{sep=',' phoenix_tsv_summaries}" ] && [ -z "~{sep=',' griphin_xlsx_summaries}" ] && [ -z "~{sep=',' griphin_tsv_summaries}" ] && [ -z "~{sep=',' ncbi_sra_excel_files}" ] && [ -z "~{sep=',' ncbi_biosample_excel_files}" ]; then
    echo "ERROR: No summary files were passed, please pick an array of files to combine."
    ls
    exit 1
  #check that something was made. If no files were created fail to let to user know
  elif [ ! -s "${combined_phoenix_tsv_summary_name}" ] && [ ! -s "${combined_griphin_tsv_summary_name}" ] && [ ! -s "${combined_griphin_xlsx_summary_name}" ] && [ ! -s "${combined_ncbi_sra_summary_name}" ] && [ ! -s "${combined_ncbi_biosample_summary_name}" ]; then
    echo "ERROR: No summary files were created something went wrong."
    ls
    exit 1
  fi

  >>>
  output {
    File?   phoenix_tsv_summary     = if defined(combined_phoenix_tsv_prefix) then "~{combined_phoenix_tsv_prefix}_Phoenix_Summary.tsv" else "Phoenix_Summary.tsv"
    File?   griphin_xlsx_summary    = if defined(combined_griphin_xlsx_prefix) then "~{combined_griphin_xlsx_prefix}_GRiPHin_Summary.xlsx" else "GRiPHin_Summary.xlsx"
    File?   griphin_tsv_summary     = if defined(combined_griphin_tsv_prefix) then "~{combined_griphin_tsv_prefix}_GRiPHin_Summary.tsv" else "GRiPHin_Summary.tsv"
    File?   biosample_excel_summary = if defined(combined_ncbi_biosample_xlsx_prefix) then "~{combined_ncbi_biosample_xlsx_prefix}_BiosampleAttributes_Microbe.1.0.xlsx" else "BiosampleAttributes_Microbe.1.0.xlsx"
    File?   sra_excel_summary       = if defined(combined_ncbi_sra_xlsx_prefix) then "~{combined_ncbi_sra_xlsx_prefix}_Sra_Microbe.1.0.xlsx" else "Sra_Microbe.1.0.xlsx"
    String  phoenix_version         = read_string("VERSION")
    String  phoenix_docker          = "quay.io/jvhagey/phoenix:2.0.2"
    String  analysis_date           = read_string("DATE")
  }
  runtime {
    docker: "quay.io/jvhagey/phoenix@sha256:d41682797fd662a4430a0f624475b0761a94611184f26f8acab769d3263b4153" # 2.2.0
    memory: "8 GB"
    cpu: 1
    disks:  "local-disk 100 SSD"
    maxRetries: 0
    preemptible: 0
  }
}