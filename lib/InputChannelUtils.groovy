// lib/InputChannelUtils.groovy
// Place this file at <project_root>/lib/InputChannelUtils.groovy
// Nextflow automatically loads all .groovy files in lib/ at startup — no import needed.

class InputChannelUtils {

    static Object previous_updater_check(meta, ar_file, db_path, type, mode_upper, Object forceUpdate = false) {
        def isList = ar_file instanceof List
        def patterns = [
            gamma           : /ResGANNCBI_(\d{8})_srst2\.fasta/,
            amrfinder       : /amrfinderdb_v\d+\.\d+_(\d{8})\.\d+\.tar\.gz/,
            ncbi_stats_ratio: /_Assembly_stats_(\d{8})\.txt/,
            ncbi_stats_gc   : /_Assembly_stats_(\d{8})\.txt/,
            srst2           : /ResGANNCBI_(\d{8})_srst2\.fasta/,
            gamma_pf        : /PF-Replicons_(\d{8})\.fasta/
        ]
        def dbFilename = db_path.toString().tokenize('/').last()
        def match = (dbFilename =~ patterns[type])
        def ardbDate = match ? match[0][1] : null

        if (ardbDate == null) {
            System.err.println "    ERROR: Could not extract date from DB filename: ${dbFilename}"
            return [meta, finalize(ar_file instanceof List ? ar_file.flatten() : [ar_file]), forceUpdate as boolean]
        }
        
        def resultFile = null
        def needsUpdate = forceUpdate as boolean 

        
        if (mode_upper == "UPDATE_PHOENIX") {
            def allFiles = (isList ? ar_file.flatten() : [ar_file]).findAll { it != null }

            // Define dated vs undated
            def datedFiles = allFiles.findAll { it.getName() =~ /\d{8}/ }
            def undatedFiles = allFiles.findAll { !(it.getName() =~ /\d{8}/) && it.getName() =~ /all_genes\.tsv$/ }

            // Combine them, but prioritize dated files for the 'newest' check
            // If it's amrfinder, we want to consider both in the pool
            def filesToProcess = (type == "amrfinder") ? (datedFiles + undatedFiles) : datedFiles

            if (!filesToProcess) {
                return [meta, null, forceUpdate as boolean]
            }

            // Identify the file matching current DB
            def newestFile = filesToProcess.find { it.getName().contains(ardbDate) }
            def hasNewest = newestFile != null
            
            // oldFiles should be everything EXCEPT the newest one
            // We sort them so that the "next newest" is at the top if needed
            def oldFiles = filesToProcess.findAll { newestFile == null || it != newestFile }
                                        .sort { a, b -> b.getName() <=> a.getName() } // Basic sort to get most recent

            // For amrfinder, if all dated files are current and an undated file exists, use it as the old file
            if (type == "amrfinder" && oldFiles.isEmpty() && !undatedFiles.isEmpty()) {
                oldFiles = undatedFiles
            }

            if (type == "gamma") {
                needsUpdate = !hasNewest
                resultFile = hasNewest ? null : finalize(oldFiles) 
            } else {
                if (needsUpdate || !hasNewest) {
                    resultFile = finalize(oldFiles)
                    needsUpdate = true 
                } else {
                    resultFile = null
                }
            }
        } else {
            resultFile = finalize(ar_file)
            needsUpdate = false
        }

        return [meta, resultFile, needsUpdate]
    }

    // THIS MUST BE A SEPARATE STATIC METHOD
    static Object finalize(files) {
        if (files instanceof List) {
            if (!files) return null
            if (files.size() == 1) return files[0]
            
            // Sort by date YYYYMMDD descending
            return files.sort { a, b ->
                def dateA = (a.name =~ /\d{8}/) ? (a.name =~ /\d{8}/)[0].toInteger() : 0
                def dateB = (b.name =~ /\d{8}/) ? (b.name =~ /\d{8}/)[0].toInteger() : 0
                return dateB <=> dateA
            }[0]
        }
        return files ?: null
    }


    static def create_samplesheet_meta(LinkedHashMap row) {
        def meta = [:]
        if (row.directory) {
            meta.project_id      = row.directory.toString().split('/')[-2]
            meta.full_project_id = new File(row.directory).getParent().replace("[", "")
        } else if (row.fastq_1) {
            meta.full_project_id = new File(row.fastq_1).parent
            meta.project_id      = new File(row.fastq_1).parentFile.name
        }
        return meta
    }

    static def transformSamplesheets(input_ch) {
        def meta = [:]
        def matcher = input_ch.name =~ /valid_(.+)\.csv$/
        meta.project_id = matcher ? matcher[0][1] : "unknown"
        return [meta, input_ch]
    }

    static def check_file_integrity(LinkedHashMap row) {
        def meta = [:]
        meta.id = row.sample
        def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
        meta.project_id = new File(clean_path).getParent()
        def regexPattern = ".*_summary\\.txt"
        File dir = new File(clean_path + "/file_integrity/")
        def files = dir.listFiles { file -> file.name ==~ /${regexPattern}/ }
        return files ? [meta, clean_path, true] : [meta, clean_path, false]
    }

    static def get_ids(dir) {
        List<String> dirNames = []
        if (dir.isDirectory()) {
            dir.eachFile { file ->
                if (file.isDirectory() && file.listFiles().any { it.name.endsWith('_summaryline.tsv') }) {
                    dirNames << file.name
                }
            }
        }
        return dirNames
    }

    static def modifiedFileChannel(input_ch, old_string, new_string) {
        def newFileName = input_ch[1].getName().replaceAll(old_string, new_string)
        def newFilePath = input_ch[1].getParent() ? input_ch[1].getParent().resolve(newFileName) : newFileName
        return [input_ch[0], newFilePath]
    }

//    static def create_summary_files_channels(LinkedHashMap row) {
//       def meta = [:]
//        meta.id = row.sample
//        def project_id = row.directory.toString().split('/')[-2]
//        def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
//        def cleaned_path = new File(clean_path).getParent()
//        meta.project_id = cleaned_path

    static def create_summary_files_channels(LinkedHashMap row) {
        def meta = [:]
        meta.id = row.sample
        
        // Get the directory path
        def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
        def cleaned_path = new File(clean_path).getParent()
        
        // FIX: Set project_id to just the folder name, not the whole path
        //meta.project_id = cleaned_path.split('/')[-1]
        meta.project_id = cleaned_path  // full path instead of split[-1]
        def project_id = row.directory.toString().split('/')[-2]
        def software_versions      = cleaned_path + "/pipeline_info/"
        def griphin_summary_tsv    = cleaned_path + "/" + project_id + "_GRiPHin_Summary.tsv"
        def griphin_summary_excel  = cleaned_path + "/" + project_id + "_GRiPHin_Summary.xlsx"
        def phx_summary            = cleaned_path + "/Phoenix_Summary.tsv"
        
        // --- NEW LOGIC FOR UPDATER ---
        def updater_path = cleaned_path + "/update_pipeline_info/software_versions.yml"
        def updater_file = new File(updater_path)
        
        // If it exists, return the path string; if not, return an empty list []
        // Nextflow treats [] as an "optional/missing" path in a tuple
        def updater_software_versions = updater_file.exists() ? updater_file.toPath() : []
        
        return [meta, griphin_summary_excel, griphin_summary_tsv, phx_summary, software_versions, updater_software_versions]
    }

    static def create_dir_channels(LinkedHashMap row) {
        def meta = [:]
        meta.id = row.sample
        def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
        def cleaned_path = new File(clean_path).getParent()
        meta.project_id = cleaned_path
        return [meta, cleaned_path]
    }

    static def add_entry_meta(input_ch) {
        def meta = [:]
        meta.project_id = input_ch[0].project_id
        def file = input_ch[1]
        meta.entry = file.text.contains('BUSCO')
        return [meta, input_ch[1]]
    }

    static def create_groups_id_and_busco(input_ch, project_folder) {
        def meta = [:]
        meta.project_id = project_folder.toString().split('/')[-1]
        meta.entry = input_ch[1].text.contains('BUSCO')
        return [meta, input_ch[1]]
    }

    static def create_sample_dir_channels(LinkedHashMap row) {
        def meta = [:]
        meta.id = row.sample
        def clean_path = row.directory.toString().endsWith("/") ? row.directory.toString()[0..-2] : row.directory.toString()
        meta.project_id = new File(clean_path).getParent()
        return [meta, clean_path]
    }

    static def create_groups_and_id(input_ch, project_folder) {
        def meta = [:]
        meta.id = input_ch[0]
        meta.project_id = project_folder.toString().split('/')[-1]
        return [meta, input_ch[1]]
    }

    static def create_groups(input_ch, project_folder) {
        def meta = [:]
        meta.id = input_ch[0]
        meta.project_id = project_folder.toString()
        return [meta, project_folder]
    }

    static def get_only_passing_samples(filePaths) {
        def passingFiles = []
        filePaths.each { filePath ->
            if (filePath.toString().contains('.nextflow')) {
                println("Skipping: ${filePath} (contains .nextflow dir).")
                return
            }
            def fileObj = new File(filePath.toString())
            if (!fileObj.text.contains('FAILED') || fileObj.text.contains('End_of_File')) {
                def fileName = filePath.getFileName().toString().split('/')[-1]
                def passing = fileName
                    .replaceAll(/_scaffolds_summary\.txt/, '')
                    .replaceAll(/_rawstats_summary\.txt/, '')
                    .replaceAll(/_corruption_summary\.txt/, '')
                    .replaceAll(/_trimstats_summary\.txt/, '')
                    .replaceAll(/_summary\.txt/, '')
                passingFiles << passing
            }
        }
        return passingFiles
    }

    static def get_old_fairy_samples(filePaths) {
        def passingFiles = []
        filePaths.each { filePath ->
            if (filePath.toString().contains('.nextflow')) {
                println("Skipping: ${filePath} (contains .nextflow dir).")
                return
            }
            def fileObj = new File(filePath.toString())
            if (fileObj.exists() && fileObj.readLines().size() >= 5) {
                def fileName = filePath.getFileName().toString().split('/')[-1]
                def passing = fileName
                    .replaceAll(/_scaffolds_summary\.txt/, '')
                    .replaceAll(/_rawstats_summary\.txt/, '')
                    .replaceAll(/_corruption_summary\.txt/, '')
                    .replaceAll(/_trimstats_summary\.txt/, '')
                    .replaceAll(/_summary\.txt/, '')
                passingFiles << passing
            }
        }
        return passingFiles
    }

    static def get_failed_samples(filePaths) {
        def failedFiles = []
        filePaths.each { filePath ->
            if (filePath.toString().contains('.nextflow')) {
                println("Skipping: ${filePath} (contains .nextflow dir).")
                return
            }
            def fileObj = new File(filePath.toString())
            if (fileObj.text.contains('FAILED') || fileObj.text.contains('End_of_File')) {
                def fileName = filePath.getFileName().toString().split('/')[-1]
                def failed = fileName
                    .replaceAll(/_scaffolds_summary\.txt/, '')
                    .replaceAll(/_rawstats_summary\.txt/, '')
                    .replaceAll(/_corruption_summary\.txt/, '')
                    .replaceAll(/_trimstats_summary\.txt/, '')
                    .replaceAll(/_summary\.txt/, '')
                failedFiles << failed
            }
        }
        return failedFiles
    }

    static def append_to_path(full_path, string) {
        return full_path.toString().endsWith('/') ? full_path.toString() + string : full_path.toString() + '/' + string
    }

    static def create_meta_with_wildcard(sample, file_extension, indir) {
        def meta = [:]
        meta.id         = sample.getName().replaceAll(file_extension, "").split('_REFSEQ_\\d{8}')[0]
        meta.project_id = indir.toString().split('/')[-1]
        return [meta, sample]
    }

    static def create_meta(sample, file_extension, indir, extra_check) {
        def meta = [:]
        if (extra_check == true) {
            def full_ext1 = "_scaffolds" + file_extension
            def full_ext2 = "_corruption" + file_extension
            def full_ext3 = "_trimstats" + file_extension
            def full_ext4 = "_rawstats" + file_extension
            meta.id = sample.getName()
                .replaceAll(full_ext1, '').replaceAll(full_ext2, '')
                .replaceAll(full_ext3, '').replaceAll(full_ext4, '')
                .replaceAll(file_extension, '')
        } else {
            meta.id = sample.getName().replaceAll(file_extension, "")
            if (file_extension == ".filtered.scaffolds.fa.txt") {
                meta.id = meta.id.tokenize('.').last()
            }
        }
        meta.project_id = indir.toString().split('/')[-1]
        return [meta, sample]
    }

    static def create_meta_non_extension(sample, indir) {
        def meta = [:]
        meta.id = sample.getSimpleName()
        def hyperVirulencePattern = /_HyperVirulence_\d{8}/
        def pfRepliconsPattern    = /_PF-Replicons_\d{8}/
        def assemblyratioPattern  = /_Assembly_ratio_\d{8}/
        def gcPattern             = /_GC_content_\d{8}/
        def arPattern             = /_ResGANNCBI_\d{8}_srst2/
        def amrfinderPattern      = /_all_genes/
        def srst2Pattern          = /__fullgenes___results/
        if (meta.id =~ hyperVirulencePattern || meta.id =~ pfRepliconsPattern ||
            meta.id =~ assemblyratioPattern  || meta.id =~ gcPattern          ||
            meta.id =~ arPattern             || meta.id =~ amrfinderPattern) {
            meta.id = meta.id
                .replaceAll(hyperVirulencePattern, '').replaceAll(pfRepliconsPattern, '')
                .replaceAll(assemblyratioPattern,  '').replaceAll(gcPattern, '')
                .replaceAll(arPattern,             '').replaceAll(amrfinderPattern, '')
                .replaceAll(srst2Pattern, '').trim()
        }
        meta.project_id = indir.toString().split('/')[-1]
        return [meta, sample]
    }

    static def check_scaffolds(scaffold_channel) {
        if (scaffold_channel[1].toString().endsWith(".fasta.gz") || scaffold_channel[1].toString().endsWith(".fa.gz")) {
            return scaffold_channel
        }
        throw new Exception("ERROR: No scaffolds found in '${scaffold_channel[1]}'. Scaffold files must end in '.fa.gz' or '.fasta.gz'.")
    }
}
