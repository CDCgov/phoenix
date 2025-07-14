version 1.0

task phoenix {
  input {
    File?    read1
    File?    read2
    File?    input_assembly
    String   samplename
    String   kraken2db = "null"
    String   pipeline = "PHOENIX"
    String   scaffold_ext = ".scaffolds.fa.gz"
    Boolean? create_ncbi_sheet = false
    Boolean? centar = false
    Int?     coverage = 30
    Int      memory = 64
    Int      cpu = 8
    Int      disk_size = 100
  }
  command <<<
    date | tee DATE
    version="v2.2.0-dev" 
    echo $version | tee VERSION

    # Debug
    export TMP_DIR=$TMPDIR
    export TMP=$TMPDIR
    env

    # Convert to lowercase using bash
    pipeline_upper=$(echo ~{pipeline} | tr '[:lower:]' '[:upper:]')

    if [ ${pipeline_upper} == "SRA" ] || [ ${pipeline_upper} == "CDC_SRA" ]; then
      # Make sample form
      echo "~{samplename}" > sample.csv
      # Run PHoeNIx
      mkdir ~{samplename}
      cd ~{samplename}
      #set input variable
      #setting --use_sra to keep ~{samplename} working
      input_file="--input_sra ../sample.csv --use_sra"
      #set scaffold as blank variable
      scaffold_ext=""
    elif [ ${pipeline_upper} == "SCAFFOLDS" ] || [ ${pipeline_upper} == "CDC_SCAFFOLDS" ]; then
      # Make sample form
      echo "sample,assembly" > sample.csv
      echo "~{samplename},~{input_assembly}" >> sample.csv
      # Run PHoeNIx
      mkdir ~{samplename}
      cd ~{samplename}
      #set input variable
      input_file="--input ../sample.csv"
      #set scaffold variable
      scaffold_ext="--scaffold_ext ~{scaffold_ext}"
    else
      # Make sample form
      echo "sample,fastq_1,fastq_2" > sample.csv
      echo "~{samplename},~{read1},~{read2}" >> sample.csv
      # Run PHoeNIx
      mkdir ~{samplename}
      cd ~{samplename}
      #set input variable
      input_file="--input ../sample.csv"
      #set scaffold as blank variable
      scaffold_ext=""
    fi

    # set shigapass db path
    shigapass_db="/opt/conda/envs/phoenix/share/shigapass-1.5.0/db"

    #checking variables
    echo $version
    echo $input_file
    echo $scaffold_ext
    echo $create_ncbi_sheet
    echo $centar

    if nextflow run cdcgov/phoenix -plugins nf-google@1.1.3 -profile terra -r $version --pipeline ${pipeline_upper} --outdir ./phx_output --terra true $input_file --kraken2db ~{kraken2db} --coverage ~{coverage} --tmpdir $TMPDIR --max_cpus ~{cpu} --max_memory '~{memory}.GB' ~{true='--centar' false='' centar} $scaffold_ext ~{true='--create_ncbi_sheet' false='' create_ncbi_sheet} --shigapass_database $shigapass_db; then
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

    save=$(find  ./ -path "*call-phoenix*" | sed 's/.*\(gs:\/\/.*\/call-phoenix\).*/\1/')
    echo $save

    find  ./ -path "*call-phoenix*" | sed 's/.*\(gs:\/\/.*\/call-phoenix\).*/\1/' | tee PROJECT_DIR

    #sed 's/.*\(gs:\/\/.*\/call-phoenix\).*/\1/' phoenix.log | sort -u | tee PROJECT_DIR

    # Get N50 from Quast file
    grep '^N50' ~{samplename}/phx_output/~{samplename}/quast/~{samplename}_summary.tsv | awk -F '\t' '{print $2}' | tee N50

    # Get AMRFinder+ output
    awk -F '\t' 'BEGIN{OFS=":"} {print $7,$12}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_AMR_CLASSES
    awk -F '\t' '{ if($8 == "core") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_AMR_CORE_GENES
    awk -F '\t' '{ if($8 == "plus") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_AMR_PLUS_GENES
    awk -F '\t' 'BEGIN{OFS=":"} {print $7,$13}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv | tail -n+2 | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_AMR_SUBCLASSES
    awk -F '\t' '{ if($9 == "STRESS") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_STRESS_GENES
    awk -F '\t' '{ if($9 == "VIRULENCE") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_VIRULENCE_GENES
    awk -F '\t' '{ if($11 == "BETA-LACTAM") { print $6}}' ~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv | tr '\n' ', ' | sed 's/.$//' | tee AMRFINDERPLUS_BETA_LACTAM_GENES

    # Gather Phoenix Output
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f5 | tee QC_OUTCOME
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f6 | tee QC_ISSUES
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f7 | awk -F',' '{print NF}' | tee WARNING_COUNT
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f7 | tee WARNINGS
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f14 | tee ESTIMATED_COVERAGE
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f17 | tee GENOME_LENGTH
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f18 | tee ASSEMBLY_RATIO
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f19 | tee ASSEMBLY_RATIO_STDEV
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f16 | tee NUM_SCAFFOLDS
    sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f15 | tee GC_PERCENT
    #sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f2,3 | tr '\t' '/' | tee PROJECT_DIR
    if [ ${pipeline_upper} == "PHOENIX" ] || [ ${pipeline_upper} == "SRA" ] || [ ${pipeline_upper} == "SCAFFOLDS" ]; then
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f20 | tee FINAL_TAXA_ID
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f21 | tee TAXA_SOURCE
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f26 | tee FASTANI_CONFIDENCE
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f27 | tee FASTANI_COVERAGE
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f25 | tee FASTANI_TAXA
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f24 | tee SHIGAPASS_TAXA
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f22 | tee KRAKEN2_TRIMD
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f23 | tee KRAKEN2_WEIGHTED
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f29 | tee MLST_SCHEME_1
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f31 | tee MLST_1
      # handling for abaumannii and ecoli primary schemes, novels
      if [[ "$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $29); print $31 "_" $29}')"=="-_" ]]; then
        sed -n 2p phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $31); print "ML" $29}' | tee MLST1_NCBI
      elif [[ "$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $29); print $31  "_" $29}')"=="*Novel*" ]]; then
        echo "" | tee MLST1_NCBI
      else
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $29); print "ML" $31 "_" $29}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST1_NCBI
      fi
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f31 | tee MLST_SCHEME_2
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f33 | tee MLST_2
      # handling for abaumannii and ecoli primary schemes, novels
      if [[ "$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $33); print $31 "_" $33}')"=="-_" ]]; then
        sed -n 7p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $33); print "ML" $31}' | tee MLST2_NCBI
      elif [[ "$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $33); print $31 "_" $33}')"=="*Novel*" ]]; then
        echo "" | tee MLST2_NCBI
      else
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $33); print "ML" $31 "_" $33}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST2_NCBI
      fi
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f22 | tee GAMMA_BETA_LACTAM_RESISTANCE_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f23 | tee GAMMA_OTHER_AR_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f24 | tee AMRFINDER_POINT_MUTATIONS
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f25 | tee HYPERVIRULENCE_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f26 | tee PLASMID_INCOMPATIBILITY_REPLICONS
      echo "Only run in CDC_PHOENIX pipeline" | tee BUSCO_DB
      echo "Only run in CDC_PHOENIX pipeline" | tee BUSCO
    elif [ ${pipeline_upper} == "CDC_PHOENIX" ] || [ ${pipeline_upper} == "CDC_SRA" ] || [ ${pipeline_upper} == "CDC_SCAFFOLDS" ]; then
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f23 | tee BUSCO
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f22 | tee BUSCO_DB
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f20 | tee FINAL_TAXA_ID
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f21 | tee TAXA_SOURCE
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f24 | tee KRAKEN2_TRIMD
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f25 | tee KRAKEN2_WEIGHTED
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f26 | tee SHIGAPASS_TAXA
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f27 | tee FASTANI_TAXA
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f28 | tee FASTANI_CONFIDENCE
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f29 | tee FASTANI_COVERAGE
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f31 | tee MLST_SCHEME_1
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f33 | tee MLST_1
      # handling for abaumannii and ecoli 2nd schemes, novels
      if [[ "$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $31); print $33 "_" $31}')"=="-_" ]]; then
        sed -n 7p phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $31); print "ML" $33}' | tee MLST1_NCBI
      elif [[ "$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $31); print $33  "_" $31}')"=="*Novel*" ]]; then
        echo "" | tee MLST1_NCBI
      else
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $31); print "ML" $33 "_" $31}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST1_NCBI
      fi
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f35 | tee MLST_SCHEME_2
      sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | cut -d$'\t' -f37 | tee MLST_2
      # handling for abaumannii and ecoli 2nd schemes, novels
      if [[ "$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $35); print $37 "_" $35}')"=="-_" ]]; then
        sed -n 7p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $35); print "ML" $37}' | tee MLST2_NCBI
      elif [[ "$(sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $35); print $37 "_" $35}')"=="*Novel*" ]]; then
        echo "" | tee MLST2_NCBI
      else
        sed -n 2p ~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv | awk -F'\t' '{gsub(/[^a-zA-Z0-9]/, "", $35); print "ML" $37 "_" $35}' | sed -E 's/_[^_]*(Achtman|Oxford|Pasteur)/_\1/' | tee MLST2_NCBI
      fi
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f24 | tee GAMMA_BETA_LACTAM_RESISTANCE_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f25 | tee GAMMA_OTHER_AR_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f26 | tee AMRFINDER_POINT_MUTATIONS
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f27 | tee HYPERVIRULENCE_GENES
      sed -n 2p ~{samplename}/phx_output/Phoenix_Summary.tsv | cut -d$'\t' -f28 | tee PLASMID_INCOMPATIBILITY_REPLICONS
    else
      echo "Pipeline not recognized. Enter one: PHOENIX, CDC_PHOENIX, SCAFFOLDS, CDC_SCAFFOLDS, SRA, or CDC_SRA."
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
    String  estimated_coverage                = read_string("ESTIMATED_COVERAGE") #make string for cases where it's "unknown"
    String  genome_length                     = read_string("GENOME_LENGTH") #make string for cases where it's "unknown"
    String  N50                               = read_string("N50")
    String  assembly_ratio                    = read_string("ASSEMBLY_RATIO")
    String  assembly_ratio_stdev              = read_string("ASSEMBLY_RATIO_STDEV")
    String  scaffold_count                    = read_string("NUM_SCAFFOLDS") #make string for cases where it's "unknown"
    String  gc_percent                        = read_string("GC_PERCENT") #make string for cases where it's "unknown"
    String  final_taxa_id                     = read_string("FINAL_TAXA_ID")
    String  taxa_source                       = read_string("TAXA_SOURCE")
    String  busco                             = read_string("BUSCO") #seems like I can't make these busco things optional so just leaving blank
    String  busco_db                          = read_string("BUSCO_DB")
    String  kraken2_trimmed                   = read_string("KRAKEN2_TRIMD")
    String  kraken2_weighted                  = read_string("KRAKEN2_WEIGHTED")
    String  shigapass_taxa                    = read_string("SHIGAPASS_TAXA")
    String  fastani_taxa                      = read_string("FASTANI_TAXA")
    String  fastani_confidence                = read_string("FASTANI_CONFIDENCE")
    String  fastani_coverage                  = read_string("FASTANI_COVERAGE") #make string for cases where it's "unknown"
    String  mlst_scheme_1                     = read_string("MLST_SCHEME_1")
    String  mlst_1                            = read_string("MLST_1")
    String  mlst1_ncbi                        = read_string("MLST1_NCBI")
    String  mlst_scheme_2                     = read_string("MLST_SCHEME_2")
    String  mlst_2                            = read_string("MLST_2")
    String  mlst2_ncbi                        = read_string("MLST2_NCBI")
    String  gamma_beta_lactam_genes           = read_string("GAMMA_BETA_LACTAM_RESISTANCE_GENES")
    String  gamma_other_ar_genes              = read_string("GAMMA_OTHER_AR_GENES")
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
    File full_results             = "~{samplename}.tar.gz"
    File griphin_excel_summary    = "~{samplename}/phx_output/phx_output_GRiPHin_Summary.xlsx"
    File griphin_tsv_summary      = "~{samplename}/phx_output/phx_output_GRiPHin_Summary.tsv"
    File phoenix_tsv_summary      = "~{samplename}/phx_output/Phoenix_Summary.tsv"
    #phoenix fastqc - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? raw_read1_html          = "~{samplename}/phx_output/~{samplename}/qc_stats/~{samplename}_1_fastqc.html"
    File? raw_read1_zip           = "~{samplename}/phx_output/~{samplename}/qc_stats/~{samplename}_1_fastqc.zip"
    File? raw_read2_html          = "~{samplename}/phx_output/~{samplename}/qc_stats/~{samplename}_2_fastqc.html"
    File? raw_read2_zip           = "~{samplename}/phx_output/~{samplename}/qc_stats/~{samplename}_2_fastqc.zip"
    #phoenix trimmed kraken/krona - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? kraken_trimd_summary    = "~{samplename}/phx_output/~{samplename}/kraken2_trimd/~{samplename}.kraken2_trimd.summary.txt"
    File? kraken_trimd_top_taxa   = "~{samplename}/phx_output/~{samplename}/kraken2_trimd/~{samplename}.kraken2_trimd.top_kraken_hit.txt"
    File? trimd_html              = "~{samplename}/phx_output/~{samplename}/kraken2_trimd/krona/~{samplename}_trimd.html"
    File? trimd_krona             = "~{samplename}/phx_output/~{samplename}/kraken2_trimd/krona/~{samplename}_trimd.krona"
    #File? classified_1            = "~{samplename}/phx_output/~{samplename}/kraken2_trimd/~{samplename}.classified_1.fasta.gz"
    #File? unclassified_1          = "~{samplename}/phx_output/~{samplename}/kraken2_trimd/~{samplename}.unclassified_1.fasta.gz"
    #File? classified_2            = "~{samplename}/phx_output/~{samplename}/kraken2_trimd/~{samplename}.classified_2.fasta.gz"
    #File? unclassified_2          = "~{samplename}/phx_output/~{samplename}/kraken2_trimd/~{samplename}.unclassified_2.fasta.gz"
    #phoenix QC - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    Array[File] file_integrity   = glob("~{samplename}/phx_output/~{samplename}/file_integrity/~{samplename}_*_summary.txt")
    File? paired_fastp_html       = "~{samplename}/phx_output/~{samplename}/fastp_trimd/~{samplename}.fastp.html"
    File? paired_fastp_json       = "~{samplename}/phx_output/~{samplename}/fastp_trimd/~{samplename}.fastp.json"
    File? single_fastp_html       = "~{samplename}/phx_output/~{samplename}/fastp_trimd/~{samplename}_singles.fastp.html"
    File? single_fastp_json       = "~{samplename}/phx_output/~{samplename}/fastp_trimd/~{samplename}_singles.fastp.json"
    File? trimmed_singles         = "~{samplename}/phx_output/~{samplename}/fastp_trimd/~{samplename}.singles.fastq.gz"
    File? trimmed_read1           = "~{samplename}/phx_output/~{samplename}/fastp_trimd/~{samplename}_1.trim.fastq.gz"
    File? trimmed_read2           = "~{samplename}/phx_output/~{samplename}/fastp_trimd/~{samplename}_2.trim.fastq.gz"
    File? trimmed_read_counts     = "~{samplename}/phx_output/~{samplename}/qc_stats/~{samplename}_trimmed_read_counts.txt"
    File? raw_read_counts         = "~{samplename}/phx_output/~{samplename}/raw_stats/~{samplename}_raw_read_counts.txt"
    File? adapter_removal_log     = "~{samplename}/phx_output/~{samplename}/qc_stats/~{samplename}.bbduk.log"
    #phoenix assembly - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? assembly_graph           = "~{samplename}/phx_output/~{samplename}/assembly/~{samplename}.assembly.gfa.gz"
    File? filtered_scaffolds_log   = "~{samplename}/phx_output/~{samplename}/assembly/~{samplename}.bbmap_filtered.log"
    File? contigs                  = "~{samplename}/phx_output/~{samplename}/assembly/~{samplename}.contigs.fa.gz"
    File? filtered_scaffolds       = "~{samplename}/phx_output/~{samplename}/assembly/~{samplename}.filtered.scaffolds.fa.gz"
    File? assembly_with_seq_names  = "~{samplename}/phx_output/~{samplename}/assembly/~{samplename}.renamed.scaffolds.fa.gz"
    File? assembly                 = "~{samplename}/phx_output/~{samplename}/assembly/~{samplename}.scaffolds.fa.gz"
    File? spades_log               = "~{samplename}/phx_output/~{samplename}/assembly/~{samplename}.spades.log"
    #phoenix wtasmbld kraken/krona
    File? kraken_wtasmbld_summary  = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld_weighted/~{samplename}.kraken2_wtasmbld.summary.txt"
    File? kraken_wtasmbld_top_taxa = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld_weighted/~{samplename}.kraken2_wtasmbld.top_kraken_hit.txt"
    File? wtasmbld_html            = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld_weighted/krona/~{samplename}_wtasmbld.html"
    File? wtasmbld_krona           = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld_weighted/krona/~{samplename}_wtasmbld.krona"
    ## phoenix asmbld kraken/krona -- only made when --pipeline CDC_PHOENIX or CDC_SCAFFOLDS is run
    File? kraken_asmbld_output     = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld/~{samplename}.kraken2_asmbld.classifiedreads.txt"
    File? kraken_asmbld_summary    = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld/~{samplename}.kraken2_asmbld.summary.txt"
    File? kraken_asmbld_top_taxa   = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld/~{samplename}.kraken2_asmbld.top_kraken_hit.txt"
    File? asmbld_html              = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld/krona/~{samplename}_asmbld.html"
    File? asmbld_krona             = "~{samplename}/phx_output/~{samplename}/kraken2_asmbld/krona/~{samplename}_asmbld.krona"
    #phoenix ani
    File? fast_ani                 = "~{samplename}/phx_output/~{samplename}/ANI/~{samplename}_REFSEQ_20250214.ani.txt"
    File? reformated_fast_ani      = "~{samplename}/phx_output/~{samplename}/ANI/~{samplename}_REFSEQ_20250214.fastANI.txt"
    File? top_20_taxa_matches      = "~{samplename}/phx_output/~{samplename}/ANI/mash_dist/~{samplename}_REFSEQ_20250214_best_MASH_hits.txt"
    File? mash_distance            = "~{samplename}/phx_output/~{samplename}/ANI/mash_dist/~{samplename}_REFSEQ_20250214.txt"
    #phoenix quast and mlst
    File? quast_summary            = "~{samplename}/phx_output/~{samplename}/quast/~{samplename}_summary.tsv"
    File? mlst_tsv                 = "~{samplename}/phx_output/~{samplename}/mlst/~{samplename}_combined.tsv"
    # cdc_phoenix busco and srst2 - optional for PHOENIX, SCAFFOLDS and SRA entries
    Array[File]? busco_generic    = glob("~{samplename}/phx_output/~{samplename}/BUSCO/short_summary.generic.*.filtered.scaffolds.fa.txt")
    Array[File]? busco_specific   = glob("~{samplename}/phx_output/~{samplename}/BUSCO/short_summary.specific.*.filtered.scaffolds.fa.txt")
    File? srst2                   = "~{samplename}/phx_output/~{samplename}/srst2/~{samplename}__fullgenes__ResGANNCBI_20250214_srst2__results.txt"
    #phoenix gamma
    File? gamma_ar_calls           = "~{samplename}/phx_output/~{samplename}/gamma_ar/~{samplename}_ResGANNCBI_20250519_srst2.gamma"
    File? blat_ar_calls            = "~{samplename}/phx_output/~{samplename}/gamma_ar/~{samplename}_ResGANNCBI_20250519_srst2.psl"
    File? gamma_hv_calls           = "~{samplename}/phx_output/~{samplename}/gamma_hv/~{samplename}_HyperVirulence_20220414.gamma"
    File? blat_hv_calls            = "~{samplename}/phx_output/~{samplename}/gamma_hv/~{samplename}_HyperVirulence_20220414.psl"
    File? gamma_pf_calls           = "~{samplename}/phx_output/~{samplename}/gamma_pf/~{samplename}_PF-Replicons_20250214.gamma"
    File? blat_pf_calls            = "~{samplename}/phx_output/~{samplename}/gamma_pf/~{samplename}_PF-Replicons_20250214.psl"
    #phoenix output
    File? assembly_ratio_file      = "~{samplename}/phx_output/~{samplename}/~{samplename}_Assembly_ratio_20250214.txt"
    File? gc_content_file          = "~{samplename}/phx_output/~{samplename}/~{samplename}_GC_content_20250214.txt"
    File  summary_line             = "~{samplename}/phx_output/~{samplename}/~{samplename}_summaryline.tsv"
    File  synopsis                 = "~{samplename}/phx_output/~{samplename}/~{samplename}.synopsis"
    File? best_taxa_id             = "~{samplename}/phx_output/~{samplename}/~{samplename}.tax"
    #phoenix amrfinder
    File? amrfinder_mutations      = "~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_mutations.tsv"
    File? amrfinder_taxa_match     = "~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_AMRFinder_Organism.csv"
    File? amrfinder_hits           = "~{samplename}/phx_output/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv"
    #species specific
    File? shigapass_summary       = "~{samplename}/phx_output/~{samplename}/ANI/~{samplename}_ShigaPass_summary.csv"
    File? centar_summary          = "~{samplename}/phx_output/~{samplename}/CENTAR/~{samplename}_centar_output.tsv"
    File? centar_ar_AA_gamma      = "~{samplename}/phx_output/~{samplename}/CENTAR/gamma_cdiff_specific_ar/~{samplename}_centar_ar_db_wt_AA_20240910.gamma"
    File? centar_ar_NT_gamma      = "~{samplename}/phx_output/~{samplename}/CENTAR/gamma_cdiff_specific_ar/~{samplename}_centar_ar_db_wt_NT_20240910.gamma"
    File? centar_tox_gamma        = "~{samplename}/phx_output/~{samplename}/CENTAR/gamma_cdiff_toxins/~{samplename}_Cdiff_toxins_srst2_20240909.gamma"
    File? centar_clade            = "~{samplename}/phx_output/~{samplename}/CENTAR/clade/~{samplename}_cdifficile_clade.tsv"
    # NCBI files - optional
    File? ncbi_biosample          = "~{samplename}/phx_output/*_BiosampleAttributes_Microbe.1.0.xlsx"
    File? ncbi_sra_metadata       = "~{samplename}/phx_output/*_Sra_Microbe.1.0.xlsx"
    #full results - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File versions_file            = "~{samplename}/phx_output/pipeline_info/software_versions.yml"
    File? multiqc_output          = "~{samplename}/phx_output/multiqc/multiqc_report.html"
  }
  runtime {
    docker: "quay.io/jvhagey/phoenix:2.2.0"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk ~{disk_size} SSD"
    maxRetries: 0
    preemptible: 0
  }
}
