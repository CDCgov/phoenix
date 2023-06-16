version 1.0

task phoenix {
  input {
    File?  read1
    File?  read2
    File?  input_assembly
    String samplename
    String kraken2db = "null"
    String entry = "PHOENIX"
    String docker = "quay.io/jvhagey/phoenix:2.0.0"
    String scaffold_ext = ".scaffolds.fa.gz"
    Int?   coverage = 30
    Int    memory = 64
    Int    cpu = 8
    Int    disk_size = 100
  }
  command <<<
    date | tee DATE
    echo $(nextflow pull cdcgov/phoenix 2>&1) | sed 's/^.*revision: //;' | tee VERSION

    # Debug
    export TMP_DIR=$TMPDIR
    export TMP=$TMPDIR
    env

    if [ ~{entry} == "SRA" ] || [ ~{entry} == "CDC_SRA" ]; then
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
    elif [ ~{entry} == "SCAFFOLDS" ] || [ ~{entry} == "CDC_SCAFFOLDS" ]; then
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

    if nextflow run cdcgov/phoenix -plugins nf-google@1.1.3 -profile terra -r v1.2.0-dev -entry ~{entry} --terra true $input_file --kraken2db ~{kraken2db} --coverage ~{coverage} --tmpdir $TMPDIR --max_cpus ~{cpu} --max_memory '~{memory}.GB' $scaffold_ext; then
      # Everything finished, pack up the results and clean up
      #tar -cf - work/ | gzip -n --best > work.tar.gz
      rm -rf .nextflow/ work/
      cd ..
      tar -cf - ~{samplename}/ | gzip -n --best > ~{samplename}.tar.gz
    else
      # Run failed
      tar -cf - work/ | gzip -n --best > work.tar.gz
      #save line for debugging specific file - just change "collated_versions.yml" to specific file name
      find  /cromwell_root/ -path "*work*" -name "*.command.err" | xargs -I {} bash -c "echo {} && cat {}"
      find  /cromwell_root/ -path "*work*" -name "*.command.out" | xargs -I {} bash -c "echo {} && cat {}"
      find  /cromwell_root/ -name "*.nextflow.log" | xargs -I {} bash -c "echo {} && cat {}"
      exit 1
    fi

    # Gather Phoenix Output
    sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f2 | tee QC_OUTCOME
    sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f3 | tee WARNING_COUNT
    sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f4 | tee ESTIMATED_COVERAGE
    sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f5 | tee GENOME_LENGTH
    sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f6 | tee ASSEMBLY_RATIO
    sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f7 | tee NUM_SCAFFOLDS
    sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f8 | tee GC_PERCENT
    if [ ~{entry} == "PHOENIX" ] || [ ~{entry} == "SRA" ] || [ ~{entry} == "SCAFFOLDS" ]; then
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f9 | tee SPECIES
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f10 | tee TAXA_CONFIDENCE
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f11 | tee TAXA_COVERAGE
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f12 | tee TAXA_SOURCE
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f13 | tee KRAKEN2_TRIMD
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f14 | tee KRAKEN2_WEIGHTED
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f15 | tee MLST_SCHEME_1
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f16 | tee MLST_1
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f17 | tee MLST_SCHEME_2
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f18 | tee MLST_2
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f19 | tee BETA_LACTAM_RESISTANCE_GENES
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f20 | tee OTHER_AR_GENES
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f21 | tee AMRFINDER_POINT_MUTATIONS
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f22 | tee HYPERVIRULENCE_GENES
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f23 | tee PLASMID_INCOMPATIBILITY_REPLICONS
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f24 | tee QC_REASON
    elif [ ~{entry} == "CDC_PHOENIX" ] || [ ~{entry} == "CDC_SRA" ] || [ ~{entry} == "CDC_SCAFFOLDS" ]; then
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f9 | tee BUSCO
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f10 | tee BUSCO_DB
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f11 | tee SPECIES
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f12 | tee TAXA_CONFIDENCE
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f13 | tee TAXA_COVERAGE
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f14 | tee TAXA_SOURCE
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f15 | tee KRAKEN2_TRIMD
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f16 | tee KRAKEN2_WEIGHTED
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f17 | tee MLST_SCHEME_1
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f18 | tee MLST_1
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f19 | tee MLST_SCHEME_2
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f20 | tee MLST_2
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f21 | tee BETA_LACTAM_RESISTANCE_GENES
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f22 | tee OTHER_AR_GENES
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f23 | tee AMRFINDER_POINT_MUTATIONS
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f24 | tee HYPERVIRULENCE_GENES
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f25 | tee PLASMID_INCOMPATIBILITY_REPLICONS
      sed -n 2p ~{samplename}/results/Phoenix_Output_Report.tsv | cut -d$'\t' -f26 | tee QC_REASON
    else
      echo "Entry point not recognized. Enter one: PHOENIX, CDC_PHOENIX, SCAFFOLDS, CDC_SCAFFOLDS, SRA, or CDC_SRA."
      exit 1
    fi
  >>>
  output {
    File?   work_files                        = "work.tar.gz"
    String  phoenix_version                   = read_string("VERSION")
    String  phoenix_docker                    = docker
    String  analysis_date                     = read_string("DATE")
    String  qc_outcome                        = read_string("QC_OUTCOME")
    String  warning_count                     = read_string("WARNING_COUNT")
    Float   estimated_coverage                = read_float("ESTIMATED_COVERAGE")
    Int     genome_length                     = read_int("GENOME_LENGTH")
    String  assembly_ratio                    = read_string("ASSEMBLY_RATIO")
    Int     scaffold_count                    = read_int("NUM_SCAFFOLDS")
    Float   gc_percent                        = read_float("GC_PERCENT")
    String? busco                             = read_string("BUSCO")
    String? busco_db                          = read_string("BUSCO_DB")
    String  species                           = read_string("SPECIES")
    String  taxa_confidence                   = read_string("TAXA_CONFIDENCE")
    Float   taxa_coverage                     = read_float("TAXA_COVERAGE")
    String  taxa_source                       = read_string("TAXA_SOURCE")
    String  kraken2_trimmed                   = read_string("KRAKEN2_TRIMD")
    String  kraken2_weighted                  = read_string("KRAKEN2_WEIGHTED")
    String  mlst_scheme_1                     = read_string("MLST_SCHEME_1")
    String  mlst_1                            = read_string("MLST_1")
    String  mlst_scheme_2                     = read_string("MLST_SCHEME_2")
    String  mlst_2                            = read_string("MLST_2")
    String  beta_lactam_resistance_genes      = read_string("BETA_LACTAM_RESISTANCE_GENES")
    String  other_ar_genes                    = read_string("OTHER_AR_GENES")
    String  amrfinder_point_mutations         = read_string("AMRFINDER_POINT_MUTATIONS")
    String  hypervirulence_genes              = read_string("HYPERVIRULENCE_GENES")
    String  plasmid_incompatibility_replicons = read_string("PLASMID_INCOMPATIBILITY_REPLICONS")
    String  qc_reason                         = read_string("QC_REASON")
    #summary files
    File full_results             = "~{samplename}.tar.gz"
    File griphin_report           = "~{samplename}/results/results_GRiPHin_Report.xlsx"
    #phoenix fastqc - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? raw_read1_html          = "~{samplename}/results/~{samplename}/qc_stats/~{samplename}_1_fastqc.html"
    File? raw_read1_zip           = "~{samplename}/results/~{samplename}/qc_stats/~{samplename}_1_fastqc.zip"
    File? raw_read2_html          = "~{samplename}/results/~{samplename}/qc_stats/~{samplename}_2_fastqc.html"
    File? raw_read2_zip           = "~{samplename}/results/~{samplename}/qc_stats/~{samplename}_2_fastqc.zip"
    #phoenix trimmed kraken/krona - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? kraken_trimd_output     = "~{samplename}/results/~{samplename}/kraken2_trimd/~{samplename}.kraken2_trimd.classifiedreads.txt"
    File? kraken_trimd_report     = "~{samplename}/results/~{samplename}/kraken2_trimd/~{samplename}.kraken2_trimd.report.txt"
    File? kraken_trimd_summary    = "~{samplename}/results/~{samplename}/kraken2_trimd/~{samplename}.trimd_summary.txt"
    File? trimd_html              = "~{samplename}/results/~{samplename}/kraken2_trimd/krona/~{samplename}_trimd.html"
    File? trimd_krona             = "~{samplename}/results/~{samplename}/kraken2_trimd/krona/~{samplename}_trimd.krona"
    File? classified_1            = "~{samplename}/results/~{samplename}/kraken2_trimd/~{samplename}.classified_1.fastq.gz"
    File? unclassified_1          = "~{samplename}/results/~{samplename}/kraken2_trimd/~{samplename}.unclassified_1.fastq.gz"
    File? classified_2            = "~{samplename}/results/~{samplename}/kraken2_trimd/~{samplename}.classified_2.fastq.gz"
    File? unclassified_2          = "~{samplename}/results/~{samplename}/kraken2_trimd/~{samplename}.unclassified_2.fastq.gz"
    #phoenix QC - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? paired_fastp_html       = "~{samplename}/results/~{samplename}/fastp_trimd/~{samplename}.fastp.html"
    File? paired_fastp_json       = "~{samplename}/results/~{samplename}/fastp_trimd/~{samplename}.fastp.json"
    File? single_fastp_html       = "~{samplename}/results/~{samplename}/fastp_trimd/~{samplename}_singles.fastp.html"
    File? single_fastp_json       = "~{samplename}/results/~{samplename}/fastp_trimd/~{samplename}_singles.fastp.json"
    File? trimmed_singles         = "~{samplename}/results/~{samplename}/fastp_trimd/~{samplename}.singles.fastq.gz"
    File? trimmed_read1           = "~{samplename}/results/~{samplename}/fastp_trimd/~{samplename}_1.trim.fastq.gz"
    File? trimmed_read2           = "~{samplename}/results/~{samplename}/fastp_trimd/~{samplename}_2.trim.fastq.gz"
    File? trimmed_read_counts     = "~{samplename}/results/~{samplename}/qc_stats/~{samplename}_trimmed_read_counts.txt"
    File? raw_read_counts         = "~{samplename}/results/~{samplename}/raw_stats/~{samplename}_raw_read_counts.txt"
    File? adapter_removal_log     = "~{samplename}/results/~{samplename}/qc_stats/~{samplename}.bbduk.log"
    #phoenix assembly - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? assembly_graph           = "~{samplename}/results/~{samplename}/assembly/~{samplename}.assembly.gfa.gz"
    File filtered_scaffolds_log   = "~{samplename}/results/~{samplename}/assembly/~{samplename}.bbmap_filtered.log"
    File? contigs                  = "~{samplename}/results/~{samplename}/assembly/~{samplename}.contigs.fa.gz"
    File filtered_scaffolds       = "~{samplename}/results/~{samplename}/assembly/~{samplename}.filtered.scaffolds.fa.gz"
    File assembly_with_seq_names  = "~{samplename}/results/~{samplename}/assembly/~{samplename}.renamed.scaffolds.fa.gz"
    File? assembly                 = "~{samplename}/results/~{samplename}/assembly/~{samplename}.scaffolds.fa.gz"
    File? spades_log               = "~{samplename}/results/~{samplename}/assembly/~{samplename}.spades.log"
    #phoenix wtasmbld kraken/krona
    File kraken_wtasmbld_output   = "~{samplename}/results/~{samplename}/kraken2_asmbld_weighted/~{samplename}.kraken2_wtasmbld.classifiedreads.txt"
    File kraken_wtasmbld_report   = "~{samplename}/results/~{samplename}/kraken2_asmbld_weighted/~{samplename}.kraken2_wtasmbld.report.txt"
    File kraken_wtasmbld_summary  = "~{samplename}/results/~{samplename}/kraken2_asmbld_weighted/~{samplename}.wtasmbld_summary.txt"
    File wtasmbld_html            = "~{samplename}/results/~{samplename}/kraken2_asmbld_weighted/krona/~{samplename}_wtasmbld.html"
    File wtasmbld_krona           = "~{samplename}/results/~{samplename}/kraken2_asmbld_weighted/krona/~{samplename}_wtasmbld.krona"
    #phoenix ani
    File fast_ani                 = "~{samplename}/results/~{samplename}/ANI/~{samplename}.ani.txt"
    File reformated_fast_ani      = "~{samplename}/results/~{samplename}/ANI/~{samplename}.fastANI.txt"
    File top_20_taxa_matches      = "~{samplename}/results/~{samplename}/ANI/mash_dist/~{samplename}_best_MASH_hits.txt"
    File mash_distance            = "~{samplename}/results/~{samplename}/ANI/mash_dist/~{samplename}.txt"
    #phoenix quast and mlst
    File quast_report             = "~{samplename}/results/~{samplename}/quast/~{samplename}_report.tsv"
    File mlst_tsv                 = "~{samplename}/results/~{samplename}/mlst/~{samplename}_combined.tsv"
    # cdc_phoenix busco and srst2 - optional for PHOENIX, SCAFFOLDS and SRA entries
    File? busco_generic           = glob("~{samplename}/results/~{samplename}/BUSCO/short_summary.generic.*.filtered.scaffolds.fa.txt")
    File? busco_specific          = glob("~{samplename}/results/~{samplename}/BUSCO/short_summary.specific.*.filtered.scaffolds.fa.txt")
    File? srst2                   = "~{samplename}/results/~{samplename}/srst2/~{samplename}__fullgenes__ResGANNCBI_20230517_srst2__results.txt"
    #phoenix gamma
    File gamma_ar_calls           = "~{samplename}/results/~{samplename}/gamma_ar/~{samplename}_ResGANNCBI_20230517_srst2.gamma"
    File blat_ar_calls            = "~{samplename}/results/~{samplename}/gamma_ar/~{samplename}_ResGANNCBI_20230517_srst2.psl"
    File gamma_hv_calls           = "~{samplename}/results/~{samplename}/gamma_hv/~{samplename}_HyperVirulence_20220414.gamma"
    File blat_hv_calls            = "~{samplename}/results/~{samplename}/gamma_hv/~{samplename}_HyperVirulence_20220414.psl"
    File gamma_pf_calls           = "~{samplename}/results/~{samplename}/gamma_pf/~{samplename}_PF-Replicons_20230504.gamma"
    File blat_pf_calls            = "~{samplename}/results/~{samplename}/gamma_pf/~{samplename}_PF-Replicons_20230504.psl"
    #phoenix output
    File assembly_ratio_file      = "~{samplename}/results/~{samplename}/~{samplename}_Assembly_ratio_20230504.txt"
    File gc_content_file          = "~{samplename}/results/~{samplename}/~{samplename}_GC_content_20230504.txt"
    File summary_line             = "~{samplename}/results/~{samplename}/~{samplename}_summaryline.tsv"
    File synopsis                 = "~{samplename}/results/~{samplename}/~{samplename}.synopsis"
    File best_taxa_id             = "~{samplename}/results/~{samplename}/~{samplename}.tax"
    #phoenix amrfinder
    File amrfinder_mutations      = "~{samplename}/results/~{samplename}/AMRFinder/~{samplename}_all_mutations.tsv"
    File? amrfinder_taxa_match    = "~{samplename}/results/~{samplename}/AMRFinder/~{samplename}_AMRFinder_Organism.csv"
    File amrfinder_hits           = "~{samplename}/results/~{samplename}/AMRFinder/~{samplename}_all_genes.tsv"
    #full results - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File versions_file            = "~{samplename}/results/pipeline_info/software_versions.yml"
    File? multiqc_report          = "~{samplename}/results/multiqc/multiqc_report.html"
  }
  runtime {
    docker: "~{docker}"
    memory: "~{memory} GB"
    cpu: cpu
    disks:  "local-disk ~{disk_size} SSD"
    maxRetries: 0
    preemptible: 0
  }
}
