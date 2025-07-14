version 1.0

import "../tasks/task_phoenix.wdl" as phoenix_nf

workflow phoenix_workflow {
  meta {
    description: "A WDL wrapper around the qc, assembly, AR gene calls components of phoenix."
  }
  input {
    File?   read1
    File?   read2
    File?   input_assembly
    String  samplename
    String  kraken2db
    String  mode
    Int?    coverage
    String? scaffold_ext
    Boolean? create_ncbi_sheet
    Boolean? centar
  }
  call phoenix_nf.phoenix {
    input:
      read1             = read1,
      read2             = read2,
      input_assembly    = input_assembly,
      samplename        = samplename,
      kraken2db         = kraken2db,
      mode              = mode,
      coverage          = coverage,
      scaffold_ext      = scaffold_ext,
      create_ncbi_sheet = create_ncbi_sheet,
      centar            = centar
  }
  output {
    #phoenix summary output values
    File?   work_files                        = phoenix.work_files
    String  project_dir                       = phoenix.project_dir
    String  phoenix_version                   = phoenix.phoenix_version
    String  phoenix_docker                    = phoenix.phoenix_docker
    String  analysis_date                     = phoenix.analysis_date
    String  qc_outcome                        = phoenix.qc_outcome
    String  warnings                          = phoenix.warnings
    String  estimated_coverage                = phoenix.estimated_coverage #make string for cases where it's "unknown"
    String  genome_length                     = phoenix.genome_length #make string for cases where it's "unknown"
    String  n50                               = phoenix.N50
    String  assembly_ratio                    = phoenix.assembly_ratio
    String  assembly_ratio_stdev              = phoenix.assembly_ratio_stdev #make string for cases where it's "unknown"
    String  scaffold_count                    = phoenix.scaffold_count #make string for cases where it's "unknown"
    String  gc_percent                        = phoenix.gc_percent #make string for cases where it's "unknown"
    String  final_taxa_id                     = phoenix.final_taxa_id
    String  taxa_source                       = phoenix.taxa_source
    String  busco                             = phoenix.busco
    String  busco_db                          = phoenix.busco_db
    String  kraken2_trimmed                   = phoenix.kraken2_trimmed
    String  kraken2_weighted                  = phoenix.kraken2_weighted
    String  shigapass_taxa                    = phoenix.shigapass_taxa
    String  fastani_taxa                      = phoenix.fastani_taxa
    String  fastani_confidence                = phoenix.fastani_confidence
    String  fastani_coverage                  = phoenix.fastani_coverage
    String  mlst_scheme_1                     = phoenix.mlst_scheme_1
    String  mlst_1                            = phoenix.mlst_1
    String  mlst1_ncbi                        = phoenix.mlst1_ncbi
    String  mlst_scheme_2                     = phoenix.mlst_scheme_2
    String  mlst_2                            = phoenix.mlst_2
    String  mlst2_ncbi                        = phoenix.mlst2_ncbi
    String  gamma_beta_lactam_genes           = phoenix.gamma_beta_lactam_genes
    String  gamma_other_ar_genes              = phoenix.gamma_other_ar_genes
    String  amrfinder_point_mutations         = phoenix.amrfinder_point_mutations
    String  amrfinder_amr_classes             = phoenix.amrfinder_amr_classes
    String  amrfinder_amr_subclasses          = phoenix.amrfinder_amr_subclasses
    String  amrfinder_core_genes              = phoenix.amrfinder_core_genes
    String  amrfinder_plus_genes              = phoenix.amrfinder_plus_genes
    String  amrfinder_stress_genes            = phoenix.amrfinder_stress_genes
    String  amrfinder_virulence_genes         = phoenix.amrfinder_virulence_genes
    String  amrfinder_beta_lactam_genes       = phoenix.amrfinder_beta_lactam_genes
    String  hypervirulence_genes              = phoenix.hypervirulence_genes
    String  plasmid_incompatibility_replicons = phoenix.plasmid_incompatibility_replicons
    String  qc_issues                         = phoenix.qc_issues
    #summary files
    File  full_results            = phoenix.full_results
    File  griphin_excel_summary   = phoenix.griphin_excel_summary
    File  griphin_tsv_summary     = phoenix.griphin_tsv_summary
    File  phoenix_tsv_summary     = phoenix.phoenix_tsv_summary
    #phoenix fastqc - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? raw_read1_html          = phoenix.raw_read1_html           # fastqc.html
    File? raw_read1_zip           = phoenix.raw_read1_zip            # fastqc.zip
    File? raw_read2_html          = phoenix.raw_read2_html           # fastqc.html
    File? raw_read2_zip           = phoenix.raw_read2_zip            # fastqc.zip
    #phoenix trimmed kraken/krona - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? kraken_trimd_summary    = phoenix.kraken_trimd_summary     # kraken2_trimd.summary.txt
    File? kraken_trimd_top_taxa   = phoenix.kraken_trimd_top_taxa    # trimd_top_taxa.txt
    File? trimd_html              = phoenix.trimd_html               # trimd.html
    File? trimd_krona             = phoenix.trimd_krona              # trimd.krona
    ## commented otu to save space, not really needed
    #File? classified_1            = phoenix.classified_1             # classified_1.fastq.gz
    #File? unclassified_1          = phoenix.unclassified_1           # unclassified_1.fastq.gz
    #File? classified_2            = phoenix.classified_2             # classified_2.fastq.gz
    #File? unclassified_2          = phoenix.unclassified_2           # unclassified_2.fastq.gz
    #phoenix QC - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    Array[File] file_integrity    = phoenix.file_integrity           # _summary.txt
    File? paired_fastp_html       = phoenix.paired_fastp_html        # fastp.html
    File? paired_fastp_json       = phoenix.paired_fastp_json        # fastp.json
    File? single_fastp_html       = phoenix.single_fastp_html        # singles.fastp.html
    File? single_fastp_json       = phoenix.single_fastp_json        # singles.fastp.json
    File? trimmed_singles         = phoenix.trimmed_singles          # singles.fastq.gz
    File? trimmed_read1           = phoenix.trimmed_read1            # read_1.trim.fastq.gz
    File? trimmed_read2           = phoenix.trimmed_read2            # read_2.trim.fastq.gz
    File? trimmed_read_counts     = phoenix.trimmed_read_counts      # trimmed_read_counts.txt
    File? raw_read_counts         = phoenix.raw_read_counts          # raw_read_counts.txt
    File? adapter_removal_log     = phoenix.adapter_removal_log      # bbduk.log
    #phoenix assembly - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File? assembly_graph          = phoenix.assembly_graph           # gfa.gz
    File? filtered_scaffolds_log  = phoenix.filtered_scaffolds_log   # bbmap_filtered.log
    File? contigs                 = phoenix.contigs                  # contigs.fa.gz
    File? filtered_scaffolds      = phoenix.filtered_scaffolds       # filtered.scaffolds.fa.gz
    File? assembly_with_seq_names = phoenix.assembly_with_seq_names  # renamed.scaffolds.fa.gz
    File? assembly                = phoenix.assembly                 # scaffolds.fa.gz
    File? spades_log              = phoenix.spades_log               # spades.log
    #phoenix wtasmbld kraken/krona
    File? kraken_wtasmbld_summary  = phoenix.kraken_wtasmbld_summary  # kraken2_wtasmbld.summary.txt
    File? kraken_wtasmbld_top_taxa = phoenix.kraken_wtasmbld_top_taxa # wtasmbld_top_taxa.txt
    File? wtasmbld_html            = phoenix.wtasmbld_html            # wtasmbld.html
    File? wtasmbld_krona           = phoenix.wtasmbld_krona           # wtasmbld.krona
    File? kraken_asmbld_output     = phoenix.kraken_asmbld_output     # kraken2_asmbld.classifiedreads.txt 
    File? kraken_asmbld_summary    = phoenix.kraken_asmbld_summary    # kraken2_asmbld.summary.txt
    File? kraken_asmbld_top_taxa   = phoenix.kraken_asmbld_top_taxa   # wtasmbld_top_taxa.txt
    File? asmbld_html              = phoenix.asmbld_html              # wtasmbld.html
    File? asmbld_krona             = phoenix.asmbld_krona             # wtasmbld.krona
    #phoenix ani
    File? fast_ani                 = phoenix.fast_ani                 # ani.txt
    File? reformated_fast_ani      = phoenix.reformated_fast_ani      # fastANI.txt
    File? top_20_taxa_matches      = phoenix.top_20_taxa_matches      # best_MASH_hits.txt 
    File? mash_distance            = phoenix.mash_distance            # .txt
    #phoenix quast and mlst
    File? quast_summary            = phoenix.quast_summary            # _report.tsv
    File? mlst_tsv                 = phoenix.mlst_tsv                 # .tsv
    # cdc_phoenix busco and srst2 - optional for PHOENIX, SCAFFOLDS and SRA entries
    Array[File]? busco_generic    = phoenix.busco_generic            # short_summary.generic.*.filtered.scaffolds.fa.txt"
    Array[File]? busco_specific   = phoenix.busco_specific           # short_summary.specific.*.filtered.scaffolds.fa.txt"
    File? srst2                   = phoenix.srst2                    # __fullgenes__ResGANNCBI_20210507_srst2__results.txt"
    #phoenix gamma
    File? gamma_ar_calls           = phoenix.gamma_ar_calls           # ResGANNCBI_20210507_srst2.gamma
    File? blat_ar_calls            = phoenix.blat_ar_calls            # ResGANNCBI_20210507_srst2.psl
    File? gamma_hv_calls           = phoenix.gamma_hv_calls           # HyperVirulence_20220414.gamma
    File? blat_hv_calls            = phoenix.blat_hv_calls            # HyperVirulence_20220414.psl
    File? gamma_pf_calls           = phoenix.gamma_pf_calls           # PF-Replicons_20220414.gamma
    File? blat_pf_calls            = phoenix.blat_pf_calls            # PF-Replicons_20220414.psl
    #phoenix output
    File? assembly_ratio_file      = phoenix.assembly_ratio_file      # Assembly_ratio_20210819.txt
    File? gc_content_file          = phoenix.gc_content_file          # GC_content_20210819.txt
    File  summary_line             = phoenix.summary_line             # summary_line.tsv
    File  synopsis                 = phoenix.synopsis                 # synopsis
    File? best_taxa_id             = phoenix.best_taxa_id             # tax
    #phoenix AMRFinder
    File? amrfinder_mutations      = phoenix.amrfinder_mutations      # all_mutations.tsv
    File? amrfinder_taxa_match     = phoenix.amrfinder_taxa_match     # AMRFinder_Organism.csv
    File? amrfinder_hits           = phoenix.amrfinder_hits           # all_genes.tsv
    #species specific
    File? shigapass_summary       = phoenix.shigapass_summary         # *_ShigaPass_summary.csv
    File? centar_summary          = phoenix.centar_summary            # *_centar_output.tsv
    File? centar_ar_AA_gamma      = phoenix.centar_ar_AA_gamma        # *_centar_ar_db_wt_AA_20240910.gamma
    File? centar_ar_NT_gamma      = phoenix.centar_ar_NT_gamma        # *_centar_ar_db_wt_NT_20240910.gamma
    File? centar_tox_gamma        = phoenix.centar_tox_gamma          # *_Cdiff_toxins_srst2_20240909.gamma
    File? centar_clade            = phoenix.centar_clade              # *_cdifficile_clade.tsv
    #File? centar_plasmid          = phoenix.centar_plasmid            # *_plasmids.tsv
    # NCBI files - optional
    File? ncbi_biosample          = phoenix.ncbi_biosample            # BiosampleAttributes_Microbe.1.0.xlsx"
    File? ncbi_sra_metadata       = phoenix.ncbi_sra_metadata         # Sra_Microbe.1.0.xlsx"
    #run files - optional for SCAFFOLDS and CDC_SCAFFOLDS entries
    File versions_file            = phoenix.versions_file             # software_versions.yml"
    File? multiqc_output          = phoenix.multiqc_output            # multiqc_report.html"
  }
}