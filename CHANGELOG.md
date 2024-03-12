# CDCgov/phoenix: Changelog

Below are the list of changes to phx since is initial release. As fixes can take multiple commits to fix the linked commit are the point at which the fix was complete. Sometimes additional changes are needed later so commits give an approximate reference for the fix. Check commits on the specific file of interest if the commit link seems off. 

## [v1.0.0](https://github.com/CDCgov/phoenix/releases/tag/v1.0.0) (10/12/2022)

🎉First official release.🎉

[Full Changelog](https://github.com/CDCgov/phoenix/compare/1.0.0-dev...v1.0.0)

## [v1.1.0](https://github.com/CDCgov/phoenix/releases/tag/v1.1.0) (03/06/2023)

[Full Changelog](https://github.com/CDCgov/phoenix/compare/v1.0.0...v1.1.0)

**Implemented Enhancements:**  
- Default branch set to main thanks @erinyoung [#84](https://github.com/CDCgov/phoenix/pull/84).  
- Added emits to allow linking of workflows to close [#42](https://github.com/CDCgov/phoenix/issues/42) [#e32132d](https://github.com/CDCgov/phoenix/commit/e32132dffc656214a9977ab6c0b22efb86f72f6f).  
- MLST output is now scanned for completeness of profiles by consolidating any allele tags to the ST column for easier scanning as well as known paralog alleles are marked for easier identification. In CDC_PHOENIX workflow ST types are consolidated, if applicable, to show concordance bewteen tools.
- Addition of 🔥🐎🐦🔥 [GRiPhin: General Report Pipeline from PHoeNIx](https://github.com/DHQP/griphin) output to `-entry CDC_PHOENIX` [#6291e9c](https://github.com/CDCgov/phoenix/commit/6291e9c6a90d28a61fb45e708536a9588a3d47a3). This was implemented to replace common report generated internally, which is why it is only in the `-entry CDC_PHOENIX`.  
- Changes to allow relative paths for kraken2 and BUSCO database to be passed rather than it requiring it to be a full path [#ecb3618](https://github.com/CDCgov/phoenix/commit/ecb3618a71e6b06a94f6282ee7220b88912d80e7) and [#d938a64](https://github.com/CDCgov/phoenix/commit/d938a6437f3e192dbe8af648c4400011fa0744e4).  
- `Phoenix_Output_Report.tsv` now has antibiotic genes and plasmid markers filtered to ensure quality [#d0fa32c](https://github.com/CDCgov/phoenix/commit/d0fa32c511a21b21366651b28dfb1539f800e262).  
   - Plasmid markers require >=60% length and >=98% identity to be reported  
   - Antibiotic Genes require >=90% length and >=98% identity to be reported  
- AMRFinder+ point mutation are now included in `Phoenix_Output_Report.tsv` under the column `AMRFinder_Point_Mutations`.
- In determine_taxID.sh, Upper taxonomy lineage now uses NCBI names and nodes files for the ability to assign nearly all possible taxonomies compared to the very limited options with the previous taxes.csv file  

**Output File Changes:**  
- Removed spaces in header of `*_all_genes.tsv` file from AMRFinder+ output and replace with underscore to allow for more friendly parsing [#fd048d1](https://github.com/CDCgov/phoenix/commit/fd048d1a54ca262617eeef32d85cd4f47650af23).  
- Fixed error causing PROKKA output to not be in Annotation folder [#d014aa0](https://github.com/CDCgov/phoenix/commit/d014aa00b27c1fa9e2d1b1151bc7f6c44d8a82b3).  
- Added headers to 2 files: `*.fastANI.txt` and `*.wtasmbld_summary.txt`.  
- Also, added headers to `phoenix_line_summary.tsv` see [wiki](https://github.com/CDCgov/phoenix/wiki/Running-PHoeNIx#sample-specific-files) for details.  
- MLST final output that includes different headers and organization was renamed to `*_combined.tsv` which includes srst2 types, if appicable, paralog tags, and any extra allele/profile tags.
- Taxonomy file now includes NCBI TaxID at each standard level. Example species line would like like this "s:287 aeruginosa"  

**Fixed Bugs:**  
- Edit to allow nf-tower to work [#b21d61f](https://github.com/CDCgov/phoenix/commit/b21d61f269212311737dffecd54664d7c8019f09)  
- Fixed pipeline failure when prokka throws error for sample names being too long (Error: ID must <= 37 chars long) [#e48e01f](https://github.com/CDCgov/phoenix/commit/e48e01fbac298541f55e949e8e8f04396aa791e8). Now sample name length doesn't matter.  
- Fixed bug where samples wouldn't end up in the `Phoenix_Output_Report.tsv` due to srst2 not finding any AR genes so the file wasn't created. Now blank file is created and remaining sample informatin is in the `Phoenix_Output_Report.tsv` [#2f52edc](https://github.com/CDCgov/phoenix/commit/2f52edc218716404b37a1e1470234e2aa32e82b3). This change only occured in `-entry CDC_PHOENIX`.  
- Fixed issue where `cp` error was thrown when relative path was given for output directory [#0c0ca55](https://github.com/CDCgov/phoenix/commit/0c0ca554861b7da28567694adc0920a6d8046d5b) and [#d938a64](https://github.com/CDCgov/phoenix/commit/d938a6437f3e192dbe8af648c4400011fa0744e4).  
- MLST PARALOGS for *Acinetobacter baumannii* are surpressed in GRiPHin report as they are...

**Database Updates:**  
- AMRFinder+ database is now static and included in the database folder [#a5d2d03](https://github.com/CDCgov/phoenix/commit/a5d2d03be4876c73b0d116d2a641c7319bf44df0). We removed the automatic updating for more control of the pipeline and lockdown to prepare for possible CLIA requirements.  
   - Version [2022-08-09.1](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.10/) currently used to be the same as the one in the curated db.  
- Curated AR gene database was updated on 2022-09-15 (yyyy-mm-dd) which includes:
   - [AMRFinderPlus database](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/)  
      - Version [2022-08-09.1](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.10/)  
   - [ARG-ANNOT](http://backup.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/2041/arg-annot-v4-aa-may2018_doc.fasta)  
      - Latest version [NT v6 July 2019](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/)  
   - [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)  
      - Includes until 2022-08-08 [commit 39f4b26](https://bitbucket.org/genomicepidemiology/resfinder_db/commits/branch/master)  

**Container Updates:**  
- MLST updated from 2.22.1 to [2.23.0](https://github.com/tseemann/mlst/releases/tag/v2.23.0).  
- BBTools updated from 38.96 to [39.01](https://sourceforge.net/projects/bbmap/).  
- AMRFinder+ was updated from 3.10.40 to [3.10.45](https://github.com/ncbi/amr/releases/tag/amrfinder_v3.10.45).  
- Scripts the utilize the phoenix_base container were updated to `quay.io/jvhagey/phoenix:base_v1.1.0` which had the python library `xlsxwriter` added to it for [`GRiPHin.py`](https://github.com/CDCgov/phoenix/blob/v1.0.1/bin/GRiPHin.py).  

## [v1.1.1](https://github.com/CDCgov/phoenix/releases/tag/v1.1.1) (03/21/2023)

[Full Changelog](https://github.com/CDCgov/phoenix/compare/v1.1.0...v1.1.1)

**Implemented Enhancements:**
- `-entry CDC_PHOENIX` workflow checks all FASTQ files for corruption and creates a list of the checked files usng the FAIry (FASTQ file Assesment of Integrity) tool [commit 1111df8](https://github.com/CDCgov/phoenix/commit/651aafe6a9459e5471ce4e4efc164587170fee62). This is a required internal QC check.  
- Expanded MLST lookup of *Citrobacter* species complex [commit 43ea24d](https://github.com/CDCgov/phoenix/commit/43ea24d0206946eb9fc90e8303fc46353e6b719b) lists the new species.  
- Increased SPAdes CPUs to 8 and memory to 16GB in `base.config`.  

**Fixed Bugs:**  
- Fix for issue [#99](https://github.com/CDCgov/phoenix/issues/99) where first gene in ar, plasmid and hypervirulence genes didn't end up in the `*_summaryline.tsv`. This same error was in `Phoenix_summary_line.py` that caused the first sample to not be include in the final report.  
- Fixed tabulation error into `*_combined.tsv` output files that in some cases would show in `GRiPHin_Report.xlsx` output as a long singular line as the MLST type.  
- Fix for issue [#91](https://github.com/CDCgov/phoenix/issues/91) where Klebsiella MLST lookup would not properly match to the correct lookup database.  
- Fixed problem where samples that didn't create scaffolds, but created contigs didn't have species printed out in `Phoenix_Output_Report.tsv` details in [commit c7f7ea5](https://github.com/CDCgov/phoenix/commit/c7f7ea5bd42a0e2010e0b15e4b4f7e9119d394a2).  
- Fixed problem in `-entry CDC_PHOENIX` where samples that didn't create scaffolds, but created contigs or samples that failed spades completely didn't have correct columns lining up in `Phoenix_Output_Report.tsv` details in [commit d17bdda](https://github.com/CDCgov/phoenix/commit/d17bdda89cf4d89aebe02a53082e5bb72c33582f).  

## [v2.0.0](https://github.com/CDCgov/phoenix/releases/tag/v2.0.0) (07/14/2023)

[Full Changelog](https://github.com/CDCgov/phoenix/compare/v1.1.1...v2.0.0)

**Implemented Enhancements:**  
- entry point for scaffolds added using either `-entry SCAFFOLDS` or `-entry CDC_SCAFFOLDS` that runs everything post SPAdes step. New input parameters `--indir` and `--scaffold_ext` added for functionality of this entry point [commit f12da60](https://github.com/CDCgov/phoenix/commit/f12da60fc4bc18499aa020ef1fb2c13d35361bb1).  
    - Supports scaffold files from shovill, spades and unicycler.
- entry point for sra added using either `-entry SRA` or `-entry CDC_SRA`. These entry points will pull samples from SRA based on what is passed to `--input_sra`, which is a file with one SRR number per line [commit a86ad3f](https://github.com/CDCgov/phoenix/commit/a86ad3fa92e287fe2be6f9631c40f9d079c5893e).  
- Check now performed on input samplesheets to confirm the same sample id, forward read and reverse read aren't used multiple times in the samplesheet [commit fd6127f](https://github.com/CDCgov/phoenix/commit/fd6127ff091d0e455a7d553415f3a5229ab6b2ec).  
- Changed many modules to `process_single` rather than `process_low` to reduce resource requirements for these steps.  
- Updates to run PHX on nf-tower with an AWS back-end. Also, updated `tower.yml` file to have working reports.  
- AMRFinder+ was updated v3.11.11 allows point mutation calling for [Burkholderia cepacia species complex](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=87882), [Burkholderia pseudomallei species complex](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=111527), ***Serratia marcescens*** and ***Staphylococcus_pseudintermedius***.
- Argument, `--coverage` added. Can be passed to increase coverage cut off that will cause sample to fail minimum qc standards (default is 30x).
- Public Kraken2 database is required rather than requesting from sharefile. For PHoeNIx >=2.0.0 you will need to download the public Standard-8 version kraken2 database **created on or after March 14th, 2023** from [Ben Langmead's github page](https://benlangmead.github.io/aws-indexes/k2). **You CANNOT use an older version of the public kraken databases** on Ben Langmead's github page. We thank @BenLangmead and @jenniferlu717 for taking the time to include an extra file in public kraken databases created after March 14th, 2023 to allow them to work in PHoeNIx!
   - For PHoeNIx <=1.1.1 you will need to download the public Standard-8 version kraken2 database created on May 17, 2021 from [Ben Langmead's github page](https://benlangmead.github.io/aws-indexes/k2). The download link is https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz.
   - The kraken database can be passed as a uncompressed folder or just in its downloaded `.tar.gz` form.

**Output File Changes:**  
- The folder `fastqc` was changed to `fastqc_trimd` to clarify it contains results from the trimmed data.  
- PROKKA module now outputs `.fsa` file (nucleotide file of genes) rather than `.fna` as the `.fna` file is really just the assembly file again.  
- Added version for base container information for `FAIRY`, `ASSET_CHECK`, `FORMAT_ANI`, `FETCH_FAILED_SUMMARIES`, `CREATE_SUMMARY_LINE`, `GATHER_SUMMARY_LINES`, and `GENERATE_PIPELINE_STATS`. This was added to `software_versions.yml`.  
- Changing the file/folder structure of some files for clarity and to make it less cluttered:
   - Folders `Annotation` and `Assembly` were changed to `annotation` and `assembly` respectively to keep continuity.  
   - Files `kraken2_asmbld/*.unclassified.fastq.gz` and `kraken2_asmbld/*.classified.fastq.gz` were changed to `kraken2_asmbld/*.unclassified.fasta.gz` and `kraken2_asmbld/*.classified.fasta.gz` as they are actually `fasta` files.  
   - `*.fastANI.txt` --> moved from `~/ANI/fastANI` to `~/ANI`.  
   - The file `*_trimmed_read_counts.txt` that was in `fastp_trimd` was moved to the folder `qc_stats`.  
   - Files `*_fastqc.zip` and `*_fastqc.html` in folder `fastqc_trimd` moved to `qc_stats`.  
   - `*.bbduk.log` --> moved from `~/removedAdapters` to `~/${sample}/qc_stats` and `removedAdapters` is not longer and output folder.  
   - `raw_stats` folder was created and contains `${sample}_raw_read_counts.txt` and `${sample}_FAIry_synopsis.txt`, previously these were in the folders `fastp_trimd` and `FAIry`, respectively.  
- Sample GC% added to `*_GC_content_20230504.txt` file.
- `*_trimmed_read_counts.txt` has `Paired_Sequenced_[reads]` column added as `Total_Sequenced_[reads]` is the number of the paired sequences and singletons.
- Files produced from FastANI, MASH and FORMAT_ANI had mash database's data appended to the file name for tracking and validation. Files are now named `*${sample}_REFSEQ_20230504.ani.txt`, `${samplename}_REFSEQ_20230504.fastANI.txt`, `${samplename}_REFSEQ_20230504_best_MASH_hits.txt` and `${samplename}_REFSEQ_20230504.txt`.
- GRiPHin file updates
   - New columns for `WARNINGS`, `ALERTS`, `Minimum_QC_Issues`, `Total_Raw_[reads]`, `Paired_Trimmed_[reads]` and `GC%`.  
   - New column `Primary_MLST_Source` as added to show if the assmebly (MLST program) or reads (SRST2) was used for MLST determination.  
   - `Auto_PassFail` and `PassFail_Reason` were changed to `Minimum_QC_Checks` and `Minimum_QC_Issues`, respectively. This was to clarifiy these are minimum requirements for QC.  
   - The column `Total_Sequenced_[bp]` was removed from the report for lack of utility.  
   - `Q30_R1_[%]`, `Q30_R2_[%]`, and `Total_Sequenced_[reads]` were relabelled as `Raw_Q30_R1_[%]`, `Raw_Q30_R2_[%]` and `Total_Trimmed_[reads]`, respectively for clarity.  

**Fixed Bugs:**  
- Added module `GET_RAW_STATS` to get raw stats, previously this was information was pulled from `FASTP_TRIMD` step, however, the input data here was post `BBDUK` which removes PhiX reads and adapters. Thus, the previous raw count was slightly off.  
- Fixed python version information not showing up for `GET_TAXA_FOR_AMRFINDER` and `GATHERING_TRIMD_READ_QC_STATS`. This was added to `software_versions.yml`.  
- Fixed issue where sample names with underscore it in caused incorrect parsing and contig number not showing up in GRiPHin reported genes [commit a0fdff5](https://github.com/CDCgov/phoenix/commit/a0fdff5536d72589535faa9bd790b8cb15f13ef7).  
- Fixed `AttributeError: 'DataFrame' object has no attribute 'map'` error that came up in GRiPhin step when your set of samples had both a macrolide and macrolide_lincosamide_streptogramin AR gene [commit 460bdbc](https://github.com/CDCgov/phoenix/commit/460bdbc05a7c01f5962289d6bff1ab6eb8de0214).  
- `Phoenix_Output_Report.tsv` was reporting %Coverage for FastANI in the `Taxa_Confidence` column rather than `%ID`. Now both are reported when FastANI is successful [commit 3b26fec](https://github.com/CDCgov/phoenix/commit/3b26fec9f15c20dfbbb0530bc1f901cb6e1119a9).  
- `GRiPHin_Report.xlsx` was switch from reported rounded numbers for coverage/similarity % to reporting the floor as reporting 100% when 99.5% is the actual number is misleading and doesn't alert the user to SNPs in genes. Now by switching to the floor 99.5% would be reported as 99% [commit 5477627](https://github.com/CDCgov/phoenix/commit/54776273e2b3ffa2231537173abb6decbccc573b).  
- Corrected GAMMA modules not printing the right version in the `software_version.yml` file [commit 5477627](https://github.com/CDCgov/phoenix/commit/54776273e2b3ffa2231537173abb6decbccc573b).  

**Database Updates:**  
- Curated AR gene database was updated on 2023-05-17 (yyyy-mm-dd) which includes:
   - [AMRFinderPlus database](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/)  
      - Version [2023-04-17.1](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/)  
   - [ARG-ANNOT](http://backup.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/2041/arg-annot-v4-aa-may2018_doc.fasta)  
      - Latest version [NT v6 July 2019](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/)  
   - [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)  
      - Bumped from `v2.0.0` to `v2.1.0` including until 2023-04-12 [commit f46d8fc](https://bitbucket.org/genomicepidemiology/resfinder_db/commits/branch/master).  
- Updated AMRFinder Database used by AMRFinder+ and GAMMA to [v2023-04-17.1](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/).  
- `SRST2_MLST` and `MLST` step now use the mlst_db which is provided in `~/phoenix/assests/databases` this is now static and no longer pulls updates from PubMLST.org. This will keep the pipeline running when PubMLST.org is down and keeps the schemes from changing if you run the same sample at different times. This was implemented to deal with PubMLST.org being down fairly often and with pipeline validation in mind.  

**Container Updates:**  
- AMRFinder+ was updated from 3.10.45 to [3.11.11](https://github.com/ncbi/amr/releases/tag/amrfinder_v3.11.11).  
- BUSCO was updated from 5.4.3 to [5.4.7](https://gitlab.com/ezlab/busco/-/blob/master/CHANGELOG).  
- MultiQC was updated from 1.11 to [1.14](https://github.com/ewels/MultiQC/releases/tag/v1.14).  
- MLST was updated from 2.22.1 to [2.23.0](https://github.com/tseemann/mlst/releases/tag/v2.23.0).


## [v2.0.1](https://github.com/CDCgov/phoenix/releases/tag/v2.0.1) (07/14/2023)

[Full Changelog](https://github.com/CDCgov/phoenix/compare/v2.0.0...v2.0.1)

**Implemented Enhancements:**  
- Updated nextflow tower scheme that describes inputs.

**Fixed Bugs:**  
- Typo fix and changed branch called in Terra task that caused Terra version to crash.

## [v2.0.2](https://github.com/CDCgov/phoenix/releases/tag/v2.0.2) (08/03/2023)

[Full Changelog](https://github.com/CDCgov/phoenix/compare/v2.0.1...v2.0.2)

**Implemented Enhancements:**  
- Added handling for -entry `SCAFFOLDS` and `CDC_SCAFFOLDS` to accept assemblies from tricylcer and flye [commit 31cb573](https://github.com/CDCgov/phoenix/commit/31cb573f1945b5bb955fb48f5f1856857f157799).  
- Added tsv version of GRiPHin_Summary.xlsx  

**Output File Changes:**  
- GRiPHin_samplesheet.csv changed to Directory_samplesheet.csv [commit b39d8d7](https://github.com/CDCgov/phoenix/commit/b39d8d706ccdd6a22de636bdd20b7cf188ae98f0)  
- In response to feedback from compliance program, "report" is being replaced by "summary" in file names to avoid confusion regarding the difference between public health results (i.e. summary) and diagnostic results (i.e. report) [commit b39d8d7](https://github.com/CDCgov/phoenix/commit/b39d8d706ccdd6a22de636bdd20b7cf188ae98f0)  
  - GRiPHin_Report.xlsx changed to GRiPHin_Summary.xlsx  
  - Phoenix_Output_Report.tsv changed to Phoenix_Summary.tsv  
  - quast/${samplename}_report.txt changed to quast/${samplename}_summary.tsv  
  - kraken2_trimd/${samplename}.trimd_summary.txt changed to kraken2_asmbld/${samplename}.kraken2_trimd.top_kraken_hit.txt  
  - kraken2_asmbld/${samplename}.asmbld_summary.txt changed to kraken2_asmbld/${samplename}.kraken2_asmbld.top_kraken_hit.txt  
  - kraken2_asmbld_weighted/${samplename}.wtasmbld_summary.txt changed to kraken2_asmbld/${samplename}.kraken2_wtasmbld.top_kraken_hit.txt  
  - kraken2_trimd/${samplename}.kraken2_trimd.report.txt changed to kraken2_trimd/${samplename}.kraken2_trimd.summary.txt  
  - kraken2_asmbld/${samplename}.kraken2_asmbld.report.txt changed to kraken2_asmbld/${samplename}.kraken2_asmbld.summary.txt  
  - kraken2_asmbld_weighted/${samplename}.kraken2_wtasmbld.report.txt changed to kraken2_asmbld_weighted/${samplename}.kraken2_wtasmbld.summary.txt  

**Fixed Bugs:**  
- For MLST when final alleles were assigned, PHX called 100% match despite 1 allele not being a match.  
- MLST step not using the custom database. A custom MLST container was added with this database included.  

**Container Updates:**  
- MLST version remains the same, but a custom database was added so that it no longer uses the database included in the software. Now hosted on quay.io.  
- Bumped up base container (v2.0.2) to have openpyxl module.  

## [v2.1.0](https://github.com/CDCgov/phoenix/releases/tag/v2.1.0) (02/11/2024)

[Full Changelog](https://github.com/CDCgov/phoenix/compare/v2.1.0...v2.0.2)

**Implemented Enhancements:**  
- Added handling for "unknown" assemblers in the scaffolds entry point so genomes can be downloaded from NCBI and run through PHoeNIx.
- For entry points CDC_PHOENIX or PHOENIX you can now use the argument `--create_ncbi_sheet` to generate partially filled out excel sheets for uploading to NCBI. You will still need to fill in some lab/sample specific information and review for accuracy, but this should speed up the process. **As a reminder, please do not submit raw sequencing data to the CDC HAI-Seq BioProject (531911) that are auto populated in these sheet unless you are a state public health laboratory, a CDC partner or have been directed to do so by DHQP. The BioProject accession IDs in these files are specifically designated for domestic HAI bacterial pathogen sequencing data, including from the Antimicrobial Resistance Laboratory Network (AR Lab Network), state public health labs, surveillance programs, and outbreaks. For inquiries about the appropriate BioProject location for your data, please contact HAISeq@cdc.gov.**    
- New Terra workflow for combining `Phoenix_Summary.tsv`, `GRiPHin_Summary.tsv` and `GRiPHin_Summary.xlsx` of multiple runs into one file. This workflow will also combine the NCBI excel sheets created when using the `--create_ncbi_sheet`.  
- `software_versions.yml` now contains versions for all custom scripts used in the pipeline to streamline its validation process and align it with CLIA requirements, ensuring smoother compliance.  
- MultiQC now contains graphs and data from BBDuk, FastP, Quast and Kraken. BUSCO is also part of MultiQC if the entry point runs it (i.e. CDC_* entries).  
- AMRFinder+ species that are screened for point mutations were updated with *Enterobacter asburiae*, *Vibrio vulfinicus* and *Vibrio parahaemolyticus*.  
- A check was added to ensure only SRR numbers are passed to -entry `CDC_SRA` and `SRA`.  
- After extensive QC cut off review addtional warnings and minimum QC cut-offs were added:
   - Minimum PASS/FAIL:
     -  %gt; 500 scaffolds
     - FAIry (file integrity check) - see Fixed Bugs section below for details.
   - Warnings:
     - 200-500 scaffolds -> high, but not enough for failure
     - Taxa Quality Checks:
        - FastANI Coverage <90% and Match <95%
        - For entries BUSCO <97% 
     - Contamination Checks: 
        - <70% of reads/weighted scaffolds assigned to top geneus hit.
        - Added weighted scaffold to kraken <30% unclassifed check (was just on reads before)
        - Added weighted scaffold to kraken only 1 genera >25% of assigned check (was just reads before)

**Output File Changes:**  
- The default outdir phx produces was changed. If the user doesn't pass `--outdir`, the default was changed from `results` to `phx_output`. This was changed in response to feedback from compliance program, to avoid confusion regarding the difference between public health results (i.e. summary) and diagnostic results (i.e. report).  
- The `phx_output/FAIry` folder will contain a `*_summaryline_failure.tsv` file for any isolate where file corruption was detected.  
- `*.tax` file had the NCBI assigned taxID added after the `:` for easy lookup.  

**Fixed Bugs:**  
- Updated `tower.yml` file to reflect file name changes in v2.0.2. This will enable nf-tower reports to properly show up. [commit e1b2b91](https://github.com/CDCgov/phoenix/commit/e1b2b912db48a55ba196f0038e5520372bb7e633)  
- `GRiPHin_Summary.xlsx` was highlighting coverage outside 40-100x despite `--coverage` setting, changes made to respect `--coverage` flag.  
- Added a fix to handle when auto select by the mlst script chooses the wrong taxonomy. PHoeNIx will force a rerun in cases where the taxonomy is known but initial mlst is run against incorrect scheme. Known instances found so far include: *E. coli* (Pasteur) being incorrectly indentified as *Aeromonas* and *E. coli* (Pasteur) being identified as *Klebsiella*. The scoring in the MLST program was updated and can now cause lower count perfect hits (e.g. 6 of 6 *Aeromonas* genes at 100%) to be scored higher than novel correct hits (e.g. 7 of 8 at 100%, 1 novel gene).  
- Corrected instance where, in some cases, an mlst scheme could not be determined that a proper out file was not created.
- Fixed issue with MLST where certain characters in filename would cause array index out of bounds error  
- Fixed issue where samples that failed SPAdes did not have `--coverage` parameter respected when generating synopsis file.  
- Fixed `-entry CDC_SCAFFOLDS` providing incorrect headers (missing `BUSCO` and `BUSCO_DB`).  
- Updated FAIry (file integrity check) to catch additional file integrity errors.  
   - FAIry detects and reports when:  
      - Corrupt fastq files that prevents the completion of gzip and zcat and generate a synopsis file when needed.  
      - If R1/R2 fastqs that do not have equal number of reads in the files.  
      - If there are no reads or scaffolds left after filtering and read trimming steps, respectively.  

**Container Updates:**  
- Containers are now called with their sha256 to streamline PHoeNIx's validation process and align it with CLIA requirements.  
- Containers updated to include developers bug fixes:  
  - fastp: v0.23.2 to v0.23.4 [bug fixes](https://github.com/OpenGene/fastp/releases).  
  - fastqc: v0.11.9 to v0.12.1 [bug fixes](https://github.com/s-andrews/FastQC/releases).  
  - kraken2: v2.1.2 to v2.1.3 which has [improvements on efficiency and bug fixes](https://github.com/s-andrews/FastQC/releases).  
  - fastani: v1.33 to v1.34 [bug fixes](https://github.com/ParBLiSS/FastANI/releases). Specifically, it fixed multi-threading output bugs. Output and interface of FastANI remains same as before.  
  - amrfinderplus: v3.11.11 to v3.11.26 which has [improvements on efficiency and bug fixes](https://github.com/ncbi/amr/releases/tag/amrfinder_v3.11.26).  
  - SRAtools v3.0.3 to 3.0.9 [updates and bug fixes](https://github.com/ncbi/sra-tools/blob/master/CHANGES.md).  
- Container for SRA entry steps `SRATOOLS_FASTERQDUMP` and `SRATOOLS_PREFETCH` was switched to a quay.io/biocontainers to address issues with the old container and ICA. [commit 68815e3](https://github.com/CDCgov/phoenix/commit/68815e3797c1944dcd0280ee658c79be90b63c0e#diff-cc23f3860dea73e90629539d540e72a7fc7cf9438de0e89eca1cc31a763c7b2b)  
- The srst2 container version stays the same, but it is now in a custom container built from [commit `73f885f55c748644412ccbaacecf12a771d0cae9`](https://github.com/CDCgov/phoenix/blob/3a270a41ebee127a3fde9b50014ce377b026987b/Dockerfiles/Dockerfile_srst2#L57C58-L57C98) as there has been a bug fix for a [rounding penalty to integer](https://github.com/katholt/srst2/commit/9eaedffb58c156e3b6c45c9273e163e2d401e792) without a new release. In addition, a fix was added to address issues related to [handling grepping  of '(' and ')'](https://github.com/CDCgov/phoenix/blob/3a270a41ebee127a3fde9b50014ce377b026987b/Dockerfiles/Dockerfile_srst2#L63). Hosting updated container on quay.io.  

**Database Updates:**  
- MLST database was pulled from PubMLST and updated on Jan 24th, 2024.  
- The Plasmid Replicons database was updated to include [an update to the Enterobacteriales.fsa database](https://bitbucket.org/genomicepidemiology/plasmidfinder_db/commits/81c11f4f2209ff12cb74b486bad4c5ede54418ad).  
- Curated AR gene database was updated on 2024-01-24 (yyyy-mm-dd) which includes:
   - [AMRFinderPlus database](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/)  
      - Version [2023-11-15.1](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/)  
   - [ARG-ANNOT](http://backup.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/2041/arg-annot-v4-aa-may2018_doc.fasta) hasn't changed since the last time the database was created and contains updates since version [NT v6 July 2019](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/)  
   - [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)  
      - Includes until 2024-01-28 [commit 97d1fe0cd0a119172037f6bdb29f8a1c7c6e6019](https://bitbucket.org/genomicepidemiology/resfinder_db/commits/branch/master)  

## [v2.1.1](https://github.com/CDCgov/phoenix/releases/tag/v2.1.1) (03/11/2024)

**Fixed Bugs:**
- Fix for issue [#130](https://github.com/CDCgov/phoenix/issues/130) Identified when an Isolate is incorrectly assigned to cronobacter scheme when it should have been ecloacae. Extension of larger scoring problem with MLST-2.23.0.
- Fixed [#142](https://github.com/CDCgov/phoenix/issues/142) where names with multiple instances of "R2" in their name couldn't be parsed properly and don't move past the corruption check step. [commit `7fc0ac3c026b7c12608be4dd1d3682675e31d0fe`](https://github.com/CDCgov/phoenix/commit/7fc0ac3c026b7c12608be4dd1d3682675e31d0fe)

**Container Updates:**  
- Containers updated to include developers bug fixes:  
  - amrfinderplus: v3.11.26 to v3.12.8 which has [changes on how AR genes are called](https://github.com/ncbi/amr/releases/tag/amrfinder_v3.12.8).  

**Database Updates:**  
- Curated AR gene database was updated on 2024-02-29 (yyyy-mm-dd) to include the new AMRFinder database:
   - [AMRFinderPlus database](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/)  
      - Version [2024-01-31.1](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/)  
   - [ARG-ANNOT](http://backup.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/2041/arg-annot-v4-aa-may2018_doc.fasta) and [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/) haven't changed since last version release.