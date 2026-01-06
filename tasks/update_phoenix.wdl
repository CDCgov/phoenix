version 1.0

task update_phoenix {
  input {
    String   samplename
    String   project_directory
    Int?     coverage = 30
    Int      memory = 64
    Int      cpu = 8
    Int      disk_size = 100
  }
  command <<<
    version="v2.2.0"
    echo $version | tee VERSION
    date | tee DATE

    #download phoenix code to get the script from
    nextflow clone cdcgov/phoenix -r $version ./$version/
    # Make sample form
    echo "sample,directory" > sample.csv
    echo "~{samplename},~{project_directory}" >> sample.csv
    # Run PHoeNIx
    mkdir ~{samplename}
    cd ~{samplename}
    #set input variable
    input_file="--input ../sample.csv"
    #set scaffold as blank variable

    # set shigapass db path
    shigapass_db="/opt/conda/envs/phoenix/share/shigapass-1.5.0/db"

    #checking variables
    echo $version
    echo $input_file

    if nextflow run cdcgov/phoenix -plugins nf-google@1.1.3 -profile terra -r $version --mode UPDATE_PHOENIX --outdir ./phx_output --terra true $input_file --coverage ~{coverage} --tmpdir $TMPDIR --max_cpus ~{cpu} --max_memory '~{memory}.GB' --shigapass_database $shigapass_db; then
      # Everything finished, pack up the results and clean up
      #tar -cf - work/ | gzip -n --best > work.tar.gz
      rm -rf .nextflow/ work/
      cd ..
      tar -cf - ~{samplename}/ | gzip -n --best > ~{samplename}.tar.gz
    else
      # Run failed
      tar -cf - work/ | gzip -n --best > work.tar.gz
      #save line for debugging specific file - just change "collated_versions.yml" to specific file name
      find  /mnt/disks/cromwell_root/~{samplename}/ -path "*work*" -name "*.command.err" | xargs -I {} bash -c "echo {} && cat {}"
      find  /mnt/disks/cromwell_root/~{samplename}/ -path "*work*" -name "*.command.out" | xargs -I {} bash -c "echo {} && cat {}"
      find  /mnt/disks/cromwell_root/~{samplename}/ -name "*.nextflow.log" | xargs -I {} bash -c "echo {} && cat {}"
      exit 1
    fi

    #check if this was a CDC run or just regular
    mode=$(head -n1 ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | grep -q "BUSCO_Lineage" && echo "CDC_PHOENIX" || echo "PHOENIX")

    # Get N50 from Quast file
    grep '^N50' ~{samplename}/phx_output/~{samplename}/quast/~{samplename}_summary.tsv | awk -F '\t' '{print $2}' | tee N50

    # Get AMRFinder+ output
    awk -F '\t' 'BEGIN{OFS=":"} {print $7,$12}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes_20250325.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_AMR_CLASSES
    awk -F '\t' '{ if($8 == "core") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes_20250325.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_AMR_CORE_GENES
    awk -F '\t' '{ if($8 == "plus") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes_20250325.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_AMR_PLUS_GENES
    awk -F '\t' 'BEGIN{OFS=":"} {print $7,$13}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes_20250325.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_AMR_SUBCLASSES
    awk -F '\t' '{ if($9 == "STRESS") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes_20250325.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_STRESS_GENES
    awk -F '\t' '{ if($9 == "VIRULENCE") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes_20250325.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_VIRULENCE_GENES
    awk -F '\t' '{ if($11 == "BETA-LACTAM") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes_20250325.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_BETA_LACTAM_GENES

    # Gather Phoenix Output
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f4 | tee QC_OUTCOME
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f5 | tee QC_ISSUES
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f6 | awk -F',' '{print NF}' | tee WARNING_COUNT
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f6 | tee WARNINGS
    #sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f2,3 | tr '\t' '/' | tee PROJECT_DIR
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f19 | tee FINAL_TAXA_ID
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f20 | tee TAXA_SOURCE
    sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f21 | tee GAMMA_BETA_LACTAM_RESISTANCE_GENES
    sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f22 | tee OTHER_AR_GENES
    sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f23 | tee AMRFINDER_POINT_MUTATIONS
    sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f24 | tee HYPERVIRULENCE_GENES
    sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f25 | tee PLASMID_INCOMPATIBILITY_REPLICONS
    if [ $mode == "PHOENIX"]; then
      if head -n 1 ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | grep -q "ShigaPass_Organism"; then
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f24 | tee SHIGAPASS_TAXA
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f29 | tee MLST_SCHEME_1
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f31 | tee MLST_1
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f33 | tee MLST_SCHEME_2
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f35 | tee MLST_2
        # handling for abaumannii and ecoli 1st schemes, novels
        MLST1_CHECK=$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $29); print $31 "_" $29}')
        if [[ "$MLST1_CHECK" == "-_" ]]; then
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($31 != "" && $31 != "-") print "ML" $31; else print ""}' | tee MLST1_NCBI
        elif [[ "$MLST1_CHECK" == *Novel* || "$MLST1_CHECK" == "Unknown_Unknown" ]]; then
          echo "" | tee MLST1_NCBI
        else
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($31 != "" && $31 != "-"); gsub(/[^a-zA-Z0-9]/, "", $29); print "ML" $31 "_" $29}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST1_NCBI
        fi
        # handling for abaumannii and ecoli 2nd schemes, novels
        MLST2_CHECK=$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $33); print $35 "_" $33}')
        if [[ "$MLST2_CHECK" == "-_" ]]; then
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($35 != "" && $35 != "-") print "ML" $35; else print ""}' | tee MLST2_NCBI
        elif [[ "$MLST2_CHECK" == *Novel* || "$MLST2_CHECK" == "Unknown_Unknown" ]]; then
          echo "" | tee MLST2_NCBI
        else
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv  | awk -F'\t' '{if ($35 != "" && $35 != "-"); {gsub(/[^a-zA-Z0-9]/, "", $33); print "ML" $35 "_" $33}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST2_NCBI
        fi
      else
        echo "" | tee SHIGAPASS_TAXA
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f28 | tee MLST_SCHEME_1
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f30 | tee MLST_1
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f32 | tee MLST_SCHEME_2
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f34 | tee MLST_2
        # handling for abaumannii and ecoli primary schemes, novels
        MLST1_CHECK=$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $28); print $30 "_" $28}')
        if [[ "$MLST1_CHECK" == "-_" ]]; then
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($30 != "" && $30 != "-") print "ML" $30; else print ""}' | tee MLST1_NCBI
        elif [[ "$MLST1_CHECK" == *Novel* || "$MLST1_CHECK" == "Unknown_Unknown" ]]; then
          echo "" | tee MLST1_NCBI
        else
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($30 != "" && $30 != "-") {gsub(/[^a-zA-Z0-9]/, "", $28); print "ML" $30 "_" $28} else print ""}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST1_NCBI
        fi
        # handling for abaumannii and ecoli primary schemes, novels
        MLST2_CHECK=$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $32); print $34 "_" $32}')
        if [[ "$MLST2_CHECK" == "-_" ]]; then
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($34 != "" && $34 != "-") print "ML" $34; else print ""}' | tee MLST2_NCBI
        elif [[ "$MLST2_CHECK" == *Novel* || "$MLST2_CHECK" == "Unknown_Unknown" ]]; then
          echo "" | tee MLST2_NCBI
        else
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($34 != "" && $34 != "-") {gsub(/[^a-zA-Z0-9]/, "", $32); print "ML" $34 "_" $32} else print ""}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST2_NCBI
        fi
      fi
    elif [ $mode == "CDC_PHOENIX" ]; then
      if head -n 1 ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | grep -q "ShigaPass_Organism"; then
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f26 | tee SHIGAPASS_TAXA
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f31 | tee MLST_SCHEME_1
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f33 | tee MLST_1
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f35 | tee MLST_SCHEME_2
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f37 | tee MLST_2
        # handling for abaumannii and ecoli 1st schemes, novels
        MLST1_CHECK=$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $31); print $33 "_" $31}')
        if [[ "$MLST1_CHECK" == "-_" ]]; then
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($33 != "" && $33 != "-") print "ML" $33; else print ""}' | tee MLST1_NCBI
        elif [[ "$MLST1_CHECK" == *Novel* || "$MLST1_CHECK" == "Unknown_Unknown" ]]; then
          echo "" | tee MLST1_NCBI
        else
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($33 != "" && $33 != "-") {gsub(/[^a-zA-Z0-9]/, "", $31); print "ML" $33 "_" $31} else print ""}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST1_NCBI
        fi
        # handling for abaumannii and ecoli 2nd schemes, novels
        MLST2_CHECK=$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $35); print $37 "_" $35}')
        if [[ "$MLST2_CHECK" == "-_" ]]; then
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($37 != "" && $37 != "-") print "ML" $37; else print ""}' | tee MLST2_NCBI
        elif [[ "$MLST2_CHECK" == *Novel* || "$MLST2_CHECK" == "Unknown_Unknown" ]]; then
          echo "" | tee MLST2_NCBI
        else
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($37 != "" && $37 != "-") {gsub(/[^a-zA-Z0-9]/, "", $35); print "ML" $37 "_" $35} else print ""}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST2_NCBI
        fi
      else
        echo "" | tee SHIGAPASS_TAXA
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f30 | tee MLST_SCHEME_1
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f32 | tee MLST_1
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f34 | tee MLST_SCHEME_2
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f36 | tee MLST_2
        # handling for abaumannii and ecoli 1st schemes, novels
        MLST1_CHECK=$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $30); print $32 "_" $30}')
        if [[ "$MLST1_CHECK" == "-_" ]]; then
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($32 != "" && $32 != "-") print "ML" $32; else print ""}' | tee MLST1_NCBI
        elif [[ "$MLST1_CHECK" == *Novel* || "$MLST1_CHECK" == "Unknown_Unknown" ]]; then
          echo "" | tee MLST1_NCBI
        else
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($32 != "" && $32 != "-") {gsub(/[^a-zA-Z0-9]/, "", $30); print "ML" $32 "_" $30} else print ""}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST1_NCBI
        fi
        # handling for abaumannii and ecoli 2nd schemes, novels
        MLST2_CHECK=$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $34); print $36 "_" $34}')
        if [[ "$MLST2_CHECK" == "-_" ]]; then
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($36 != "" && $36 != "-") print "ML" $36; else print ""}' | tee MLST2_NCBI
        elif [[ "$MLST2_CHECK" == *Novel* || "$MLST2_CHECK" == "Unknown_Unknown" ]]; then
          echo "" | tee MLST2_NCBI
        else
          sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{if ($36 != "" && $36 != "-") {gsub(/[^a-zA-Z0-9]/, "", $34); print "ML" $36 "_" $34} else print ""}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST2_NCBI
        fi
      fi
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f21 | tee GAMMA_BETA_LACTAM_RESISTANCE_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f22 | tee OTHER_AR_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f23 | tee AMRFINDER_POINT_MUTATIONS
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f24 | tee HYPERVIRULENCE_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f25 | tee PLASMID_INCOMPATIBILITY_REPLICONS
      #sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f26 | tee QC_REASON
    else
      echo "Pipeline mode not recognized. Enter UPDATE_PHOENIX or update_phoenix."
      exit 1   
    fi
  >>>
  output {
    File?   work_files                        = "work.tar.gz"
    String  project_dir                       = read_string("PROJECT_DIR")
    String  phoenix_version                   = read_string("VERSION")
    String  phoenix_docker                    = "quay.io/jvhagey/phoenix:2.2.0"
    String  analysis_date                     = read_string("DATE")
    String  qc_outcome                        = read_string("QC_OUTCOME")
    String  warnings                          = read_string("WARNINGS")
    String  final_taxa_id                     = read_string("FINAL_TAXA_ID")
    String  taxa_source                       = read_string("TAXA_SOURCE")
    String? shigapass_taxa                    = read_string("SHIGAPASS_TAXA")
    String  mlst_scheme_1                     = read_string("MLST_SCHEME_1")
    String  mlst_1                            = read_string("MLST_1")
    String  mlst1_ncbi                        = read_string("MLST1_NCBI")
    String  mlst_scheme_2                     = read_string("MLST_SCHEME_2")
    String  mlst_2                            = read_string("MLST_2")
    String  mlst2_ncbi                        = read_string("MLST2_NCBI")
    String  gamma_beta_lactam_genes           = read_string("GAMMA_BETA_LACTAM_RESISTANCE_GENES")
    String  other_ar_genes                    = read_string("OTHER_AR_GENES")
    String  amrfinder_point_mutations         = read_string("AMRFINDER_POINT_MUTATIONS")
    String  amrfinder_amr_classes             = read_string("AMRFINDERPLUS_AMR_CLASSES")
    String  amrfinder_amr_subclasses          = read_string("AMRFINDERPLUS_AMR_SUBCLASSES")
    String  amrfinder_core_genes              = read_string("AMRFINDERPLUS_AMR_CORE_GENES")
    String  amrfinder_plus_genes              = read_string("AMRFINDERPLUS_AMR_PLUS_GENES")
    String  amrfinder_stress_genes            = read_string("AMRFINDERPLUS_STRESS_GENES")
    String  amrfinder_virulence_genes         = read_string("AMRFINDERPLUS_VIRULENCE_GENES")
    String  amrfinder_beta_lactam_genes       = read_string("AMRFINDERPLUS_BETA_LACTAM_GENES")
    String  hypervirulence_genes              = read_string("HYPERVIRULENCE_GENES")
    String  plasmid_incompatibility_replicons = read_string("PLASMID_INCOMPATIBILITY_REPLICONS")
    String  qc_issues                         = read_string("QC_ISSUES")
    #summary files
    File updater_log              =  "phx_output/~{samplename}/~{samplename}_updater_log.tsv"
    File full_results             = "~{samplename}.tar.gz"
    File griphin_excel_summary    = "~{samplename}/phx_output/phx_output_GRiPHin_Summary.xlsx"
    File griphin_tsv_summary      = "~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv"
    File phoenix_tsv_summary      = "~{samplename}/phx_output/Phoenix_Summary.tsv"
    #phoenix ani
    File? reformated_fast_ani      = "~{samplename}/phx_output/~{samplename}/ANI/~{samplename}_REFSEQ_20250214.fastANI.txt"
    #phoenix quast and mlst
    File  mlst_tsv                 = "~{samplename}/phx_output/~{samplename}/mlst/~{samplename}_combined.tsv"
    # cdc_phoenix busco and srst2 - optional for PHOENIX, SCAFFOLDS and SRA entries
    File? srst2                    = "~{samplename}/phx_output/~{samplename}/srst2/~{samplename}__fullgenes__ResGANNCBI_20250519_srst2__results.txt"
    #phoenix gamma
    File  gamma_ar_calls           = "~{samplename}/phx_output/~{samplename}/gamma_ar/~{samplename}_ResGANNCBI_20250519_srst2.gamma"
    File  blat_ar_calls            = "~{samplename}/phx_output/~{samplename}/gamma_ar/~{samplename}_ResGANNCBI_20250519_srst2.psl"
    #phoenix output
    File  summary_line             = "~{samplename}/phx_output/~{samplename}/~{samplename}_summaryline.tsv"
    File  synopsis                 = "~{samplename}/phx_output/~{samplename}/~{samplename}.synopsis"
    File? best_taxa_id             = "~{samplename}/phx_output/~{samplename}/~{samplename}.tax"
    #phoenix amrfinder
    File  amrfinder_mutations      = "~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_mutations_20250325.tsv"
    File? amrfinder_taxa_match     = "~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_AMRFinder_Organism.csv"
    File? amrfinder_hits           = "~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes_20250325.tsv"
    #species specific
    File? shigapass_summary        = "~{samplename}/phx_output/~{samplename}/ANI/~{samplename}_ShigaPass_summary.csv"
    #full results - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File  versions_file            = "~{samplename}/phx_output/pipeline_info/software_versions.yml"
  }
  runtime {
    docker: "quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f" # base_v2.2.0
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk ~{disk_size} SSD"
    maxRetries: 0
    preemptible: 0
  }
}