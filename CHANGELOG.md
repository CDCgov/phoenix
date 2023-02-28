# CDCgov/phoenix: Changelog

Below are the list of changes to phx since is initial release. As fixes can take multiple commits to fix the linked commit are the point at which the fix was complete. Sometimes additional changes are needed later so commits give an approximate reference for the fix. Check commits on the specific file of interest if the commit link seems off. 

## [v1.0.0](https://github.com/CDCgov/phoenix/releases/tag/v1.0.0) (10/12/2022)

🎉First official release.🎉

[Full Changelog](https://github.com/CDCgov/phoenix/compare/1.0.0-dev...v1.0.0)

## [v1.1.0](add link) (02/XX/2023)

[Full Changelog](https://github.com/CDCgov/phoenix/compare/main...v1.0.1)

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
- MLST final output that includes different headers and organization was renamed to *_combined.tsv which includes srst2 types, if appicable, paralog tags, and any extra allele/profile tags

**Fixed Bugs:**  
- Edit to allow nf-tower to work [#b21d61f](https://github.com/CDCgov/phoenix/commit/b21d61f269212311737dffecd54664d7c8019f09)  
- Fixed pipeline failure when prokka throws error for sample names being too long (Error: ID must <= 37 chars long) [#e48e01f](https://github.com/CDCgov/phoenix/commit/e48e01fbac298541f55e949e8e8f04396aa791e8). Now sample name length doesn't matter.  
- Fixed bug where samples wouldn't end up in the `Phoenix_Output_Report.tsv` due to srst2 not finding any AR genes so the file wasn't created. Now blank file is created and remaining sample informatin is in the `Phoenix_Output_Report.tsv` [#2f52edc](https://github.com/CDCgov/phoenix/commit/2f52edc218716404b37a1e1470234e2aa32e82b3). This change only occured in `-entry CDC_PHOENIX`.  
- Fixed issue where `cp` error was thrown when relative path was given for output directory [#0c0ca55](https://github.com/CDCgov/phoenix/commit/0c0ca554861b7da28567694adc0920a6d8046d5b) and [#d938a64](https://github.com/CDCgov/phoenix/commit/d938a6437f3e192dbe8af648c4400011fa0744e4).  

**Database Updates:**  
- AMRFinder+ database is now static and included in the database folder [#a5d2d03](https://github.com/CDCgov/phoenix/commit/a5d2d03be4876c73b0d116d2a641c7319bf44df0). We removed the automatic updating for more control of the pipeline and lockdown to prepare for possible CLIA requirements. 

**Container Updates:**  
- MLST updated from 2.22.1 to [2.23.0](https://github.com/tseemann/mlst/releases/tag/v2.23.0).  
- BBTools updated from 38.96 to [39.01](https://sourceforge.net/projects/bbmap/).  
- AMRFinder+ was updated from 3.10.40 to [3.11.2](https://github.com/ncbi/amr/releases/tag/amrfinder_v3.11.2). Necessary because new amrfinder database could not be downloaded without it >v3.11.
- For folks running PHoeNIx on Terra the terra container had the following additions
    - anaconda::urllib2
- Scripts the utilize the phoenix_base container were updated to `quay.io/jvhagey/phoenix:base_v1.1.0` which had the python library `xlsxwriter` added to it for [`GRiPHin.py`](https://github.com/CDCgov/phoenix/blob/v1.0.1/bin/GRiPHin.py).  
