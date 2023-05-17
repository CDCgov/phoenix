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

**Output File Changes:**  
- Removed spaces in header of `*_all_genes.tsv` file from AMRFinder+ output and replace with underscore to allow for more friendly parsing [#fd048d1](https://github.com/CDCgov/phoenix/commit/fd048d1a54ca262617eeef32d85cd4f47650af23).  
- Fixed error causing PROKKA output to not be in Annotation folder [#d014aa0](https://github.com/CDCgov/phoenix/commit/d014aa00b27c1fa9e2d1b1151bc7f6c44d8a82b3).  
- Added headers to 2 files: `*.fastANI.txt` and `*.wtasmbld_summary.txt`.  
- Also, added headers to `phoenix_line_summary.tsv` see [wiki](https://github.com/CDCgov/phoenix/wiki/Running-PHoeNIx#sample-specific-files) for details.  
- MLST final output that includes different headers and organization was renamed to `*_combined.tsv` which includes srst2 types, if appicable, paralog tags, and any extra allele/profile tags.  

**Fixed Bugs:**  
- Edit to allow nf-tower to work [#b21d61f](https://github.com/CDCgov/phoenix/commit/b21d61f269212311737dffecd54664d7c8019f09)  
- Fixed pipeline failure when prokka throws error for sample names being too long (Error: ID must <= 37 chars long) [#e48e01f](https://github.com/CDCgov/phoenix/commit/e48e01fbac298541f55e949e8e8f04396aa791e8). Now sample name length doesn't matter.  
- Fixed bug where samples wouldn't end up in the `Phoenix_Output_Report.tsv` due to srst2 not finding any AR genes so the file wasn't created. Now blank file is created and remaining sample informatin is in the `Phoenix_Output_Report.tsv` [#2f52edc](https://github.com/CDCgov/phoenix/commit/2f52edc218716404b37a1e1470234e2aa32e82b3). This change only occured in `-entry CDC_PHOENIX`.  
- Fixed issue where `cp` error was thrown when relative path was given for output directory [#0c0ca55](https://github.com/CDCgov/phoenix/commit/0c0ca554861b7da28567694adc0920a6d8046d5b) and [#d938a64](https://github.com/CDCgov/phoenix/commit/d938a6437f3e192dbe8af648c4400011fa0744e4).  

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

## [v1.2.0](https://github.com/CDCgov/phoenix/releases/tag/v1.2.0) (0X/XX/2023)

[Full Changelog](https://github.com/CDCgov/phoenix/compare/v1.1.1...v1.2.0)

**Implemented Enhancements:**  
- Check performed on GET_RAW_STATS combined reads file to confirm read pair counts are identical before proceeding to BBDUK. Read pairs that do not match will not run through the remainder of the pipeline.
- entry point for scaffolds added using either `-entry SCAFFOLDS` or `-entry CDC_SCAFFOLDS` that runs everything post SPAdes step. New input parameters `--indir` and `--scaffold_ext` added for functionality of this entry point [commit f12da60](https://github.com/CDCgov/phoenix/commit/f12da60fc4bc18499aa020ef1fb2c13d35361bb1).  
    - Supports scaffold files from shovill, spades and unicycler.
- entry point for sra added using either `-entry SRA` or `-entry CDC_SRA`. These entry points will pull samples from SRA based on what is passed to `--input_sra`, which is a file with one SRR number per line [commit a86ad3f](https://github.com/CDCgov/phoenix/commit/a86ad3fa92e287fe2be6f9631c40f9d079c5893e).  
- Check now performed on input samplesheets to confirm the same sample id, forward read and reverse read aren't used multiple times in the samplesheet [commit fd6127f](https://github.com/CDCgov/phoenix/commit/fd6127ff091d0e455a7d553415f3a5229ab6b2ec).  
- Changed many modules to `process_single` rather than `process_low` to reduce resource requirements for these steps.  
- Updates to run PHX on nf-tower with an AWS back-end. Also, updated `tower.yml` file to have working reports.  
- AMRFinder+ was updated v3.11.11 allows point mutation calling for [Burkholderia cepacia species complex](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=87882), [Burkholderia pseudomallei species complex](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=111527), ***Serratia marcescens*** and ***Staphylococcus_pseudintermedius***.

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
   - `*.bbduk.log` --> moved from `~/removedAdapters` to `~/${meta.id}/qc_stats` and `removedAdapters` is not longer and output folder.  

**Fixed Bugs:**  
- Added module `GET_RAW_STATS` to get raw stats, previously this was information was pulled from `FASTP_TRIMD` step, however, the input data here was post `BBDUK` which removes PhiX reads and adapters. Thus, the previous raw count was slightly off.  
- Fixed python version information not showing up for `GET_TAXA_FOR_AMRFINDER` and `GATHERING_TRIMD_READ_QC_STATS`. This was added to `software_versions.yml`.  
- Fixed issue where sample names with underscore it in caused incorrect parsing and contig number not showing up in GRiPHin reported genes [commit a0fdff5](https://github.com/CDCgov/phoenix/commit/a0fdff5536d72589535faa9bd790b8cb15f13ef7).  
- Fixed `AttributeError: 'DataFrame' object has no attribute 'map'` error that came up in GRiPhin step when your set of samples had both a macrolide and macrolide_lincosamide_streptogramin AR gene [commit 460bdbc](https://github.com/CDCgov/phoenix/commit/460bdbc05a7c01f5962289d6bff1ab6eb8de0214).  
- `Phoenix_Output_Report.tsv` was reporting %Coverage for FastANI in the `Taxa_Confidence` column rather than `%ID`. Now both are reported when FastANI is successful [commit 3b26fec](https://github.com/CDCgov/phoenix/commit/3b26fec9f15c20dfbbb0530bc1f901cb6e1119a9).  
- `GRiPHin_Report.xlsx` was switch from reported rounded numbers for coverage/%similarity to reporting the floor as reporting 100% when 99.5% is the actual number is misleading and doesn't alert the user to SNPs in genes. Now by switching to the floor 99.5% would be reported as 99% [commit 5477627](https://github.com/CDCgov/phoenix/commit/54776273e2b3ffa2231537173abb6decbccc573b).  
- Corrected GAMMA modules not printing the right version in the `software_version.yml` file [commit 5477627](https://github.com/CDCgov/phoenix/commit/54776273e2b3ffa2231537173abb6decbccc573b). 

**Database Updates:**  
- Curated AR gene database was updated on 2023-XX-XX (yyyy-mm-dd) which includes:
   - [AMRFinderPlus database](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/)  
      - Version [2023-04-17.1](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/)  
   - [ARG-ANNOT](http://backup.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/2041/arg-annot-v4-aa-may2018_doc.fasta)  
      - Latest version [NT v6 July 2019](https://www.mediterranee-infection.com/acces-ressources/base-de-donnees/arg-annot-2/)  
   - [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder_db/src/master/)  
      - Bumped from `v2.0.0` to `v2.1.0` including until 2023-04-12 [commit f46d8fc](https://bitbucket.org/genomicepidemiology/resfinder_db/commits/branch/master).  
- Updated AMRFinder Database used by AMRFinder+ to [v2023-04-17.1](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/).  
- SRST2_MLST and MLST step now use the mlst_db which is provided in `~/phoenix/assests/databases` these are now static and nolonger pull updates from PubMLST.org. This will keep the pipeline running when PubMLST.org is down and keeps the schemes from changing if you run the same sample at different times. This was emplimented to deal with PubMLST.org being down fairly often and with pipeline validation in mind.

**Container Updates:**  
- AMRFinder+ was updated from 3.10.45 to [3.11.11](https://github.com/ncbi/amr/releases/tag/amrfinder_v3.11.11).  