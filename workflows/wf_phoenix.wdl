version 1.0

import "../tasks/task_phoenix.wdl" as phoenix_nf

workflow phoenix_workflow {
  meta {
    description: "A WDL wrapper around the qc, assembly, AR gene calls components of phoenix."
  }
  input {
    File   read1
    File   read2
    String samplename
    String kraken2db
  }
  call phoenix_nf.phoenix {
    input:
      read1          = read1,
      read2          = read2,
      samplename     = samplename,
      kraken2db      = kraken2db
  }
  output {
    #phoenix summary output values
    File?  work_files                        = phoenix.work_files
    String phoenix_version                   = phoenix.phoenix_version
    String phoenix_docker                    = phoenix.phoenix_docker
    String analysis_date                     = phoenix.analysis_date
    String qc_outcome                        = phoenix.qc_outcome
    String warning_count                     = phoenix.warning_count
    Float  estimated_coverage                = phoenix.estimated_coverage
    Int    genome_length                     = phoenix.genome_length
    String assembly_ratio                    = phoenix.assembly_ratio
    Int    scaffold_count                    = phoenix.scaffold_count
    Float  gc_percent                        = phoenix.gc_percent
    String species                           = phoenix.species
    String taxa_confidence                   = phoenix.taxa_confidence
    String taxa_source                       = phoenix.taxa_source
    String kraken2_trimmed                   = phoenix.kraken2_trimmed
    String kraken2_weighted                  = phoenix.kraken2_weighted
    String mlst_scheme_1                     = phoenix.mlst_scheme_1
    String mlst_1                            = phoenix.mlst_1
    String mlst_scheme_2                     = phoenix.mlst_scheme_2
    String mlst_2                            = phoenix.mlst_2
    String beta_lactam_resistance_genes      = phoenix.beta_lactam_resistance_genes
    String other_ar_genes                    = phoenix.other_ar_genes
    String amrfinder_point_mutations         = phoenix.amrfinder_point_mutations
    String hypervirulence_genes              = phoenix.hypervirulence_genes
    String plasmid_incompatibility_replicons = phoenix.plasmid_incompatibility_replicons
    String qc_reason                         = phoenix.qc_reason
    File full_results                        = phoenix.full_results
    #phoenix fastqc
    File raw_read1_html           = phoenix.raw_read1_html           # fastqc.html
    File raw_read1_zip            = phoenix.raw_read1_zip            # fastqc.zip
    File raw_read2_html           = phoenix.raw_read2_html           # fastqc.html
    File raw_read2_zip            = phoenix.raw_read2_zip            # fastqc.zip
    #phoenix trimmed kraken/krona
    File kraken_trimd_output      = phoenix.kraken_trimd_output      # kraken2_trimd.classifiedreads.txt 
    File kraken_trimd_report      = phoenix.kraken_trimd_report      # kraken2_trimd.report.txt
    File kraken_trimd_summary     = phoenix.kraken_trimd_summary     # trimd_summary.txt
    File trimd_html               = phoenix.trimd_html               # trimd.html
    File trimd_krona              = phoenix.trimd_krona              # trimd.krona
    File classified_1             = phoenix.classified_1             # classified_1.fastq.gz
    File unclassified_1           = phoenix.unclassified_1           # unclassified_1.fastq.gz
    File classified_2             = phoenix.classified_2             # classified_2.fastq.gz
    File unclassified_2           = phoenix.unclassified_2           # unclassified_2.fastq.gz
    #phoenix QC
    File paired_fastp_html        = phoenix.paired_fastp_html        # fastp.html
    File paired_fastp_json        = phoenix.paired_fastp_json        # fastp.json
    File single_fastp_html        = phoenix.single_fastp_html        # singles.fastp.html
    File single_fastp_json        = phoenix.single_fastp_json        # singles.fastp.json
    File trimmed_singles          = phoenix.trimmed_singles          # singles.fastq.gz
    File trimmed_read1            = phoenix.trimmed_read1            # read_1.trim.fastq.gz
    File trimmed_read2            = phoenix.trimmed_read2            # read_2.trim.fastq.gz
    File trimmed_read_counts      = phoenix.trimmed_read_counts      # trimmed_read_counts.txt
    File raw_read_counts          = phoenix.raw_read_counts          # raw_read_counts.txt
    File adapter_removal_log      = phoenix.adapter_removal_log      # bbduk.log
    #phoenix assembly
    File assembly_graph           = phoenix.assembly_graph           # gfa.gz
    File filtered_scaffolds_log   = phoenix.filtered_scaffolds_log   # bbmap_filtered.log
    File contigs                  = phoenix.contigs                  # contigs.fa.gz
    File filtered_scaffolds       = phoenix.filtered_scaffolds       # filtered.scaffolds.fa.gz
    File assembly_with_seq_names  = phoenix.assembly_with_seq_names  # renamed.scaffolds.fa.gz
    File assembly                 = phoenix.assembly                 # scaffolds.fa.gz
    File spades_log               = phoenix.spades_log               # spades.log
    #phoenix wtasmbld kraken/krona
    File kraken_wtasmbld_output   = phoenix.kraken_wtasmbld_output   # kraken2_wtasmbld.classifiedreads.txt 
    File kraken_wtasmbld_report   = phoenix.kraken_wtasmbld_report   # kraken2_wtasmbld.report.txt
    File kraken_wtasmbld_summary  = phoenix.kraken_wtasmbld_summary  # wtasmbld_summary.txt
    File wtasmbld_html            = phoenix.wtasmbld_html            # wtasmbld.html
    File wtasmbld_krona           = phoenix.wtasmbld_krona           # wtasmbld.krona
    #phoenix ani
    File fast_ani                 = phoenix.fast_ani                # ani.txt
    File reformated_fast_ani      = phoenix.reformated_fast_ani     # fastANI.txt
    File top_20_taxa_matches      = phoenix.top_20_taxa_matches     # best_MASH_hits.txt 
    File mash_distance            = phoenix.mash_distance           # .txt
    #phoenix quast and mlst
    File quast_report             = phoenix.quast_report            # _report.tsv
    File mlst_tsv                 = phoenix.mlst_tsv                # .tsv
    #phoenix gamma
    File gamma_ar_calls           = phoenix.gamma_ar_calls         # ResGANNCBI_20210507_srst2.gamma
    File blat_ar_calls            = phoenix.blat_ar_calls          # ResGANNCBI_20210507_srst2.psl
    File gamma_hv_calls           = phoenix.gamma_hv_calls         # HyperVirulence_20220414.gamma
    File blat_hv_calls            = phoenix.blat_hv_calls          # HyperVirulence_20220414.psl
    File gamma_pf_calls           = phoenix.gamma_pf_calls         # PF-Replicons_20220414.gamma
    File blat_pf_calls            = phoenix.blat_pf_calls          # PF-Replicons_20220414.psl
    #phoenix output
    File assembly_ratio_file      = phoenix.assembly_ratio_file    # Assembly_ratio_20210819.txt
    File gc_content_file          = phoenix.gc_content_file        # GC_content_20210819.txt
    File summary_line             = phoenix.summary_line           # summary_line.tsv
    File synopsis                 = phoenix.synopsis               # synopsis
    File best_taxa_id             = phoenix.best_taxa_id           # tax
    #phoenix AMRFinder
    File amrfinder_mutations      = phoenix.amrfinder_organism      # all_mutations.tsv
    File? amrfinder_taxa_match    = phoenix.amr_taxa_match          # AMRFinder_Organism.csv
    File amrfinder_hits           = phoenix.amr_hits                # all_genes.tsv
    #run files
    File versions_file            = phoenix.versions_file           # software_versions.yml"
    File multiqc_report           = phoenix.multiqc_report          # multiqc_report.html"
  }
}