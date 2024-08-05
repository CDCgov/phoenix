//
// workflow handles taking in either a samplesheet or directory and creates correct channels for scaffolds entry point
//

include { SAMPLESHEET_CHECK  } from '../../modules/local/samplesheet_check'
include { CREATE_SAMPLESHEET } from '../../modules/local/create_samplesheet'

workflow CREATE_UPDATE_INPUT_CHANNEL {
    take:
        indir        // params.indir
        samplesheet  // params.input
        ch_versions

    main:
        //if input directory is passed use it to gather assemblies otherwise use samplesheet
        if (indir != null) {

            // get file_integrity file for MLST updating
            def scaffolds_integrity_glob = append_to_path(params.indir.toString(),'*/file_integrity/*_scaffolds_summary.txt')
            //create file_integrity file channel with meta information -- need to pass to DO_MLST subworkflow
            file_integrity_ch = Channel.fromPath(scaffolds_integrity_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_scaffolds_summary.txt")}

            //get list of all samples in the folder - just using the file_integrity file, but it could be any file really
            def file_integrity_glob = append_to_path(params.indir.toString(),'*/file_integrity/*_*_summary.txt')
            //create file_integrity file channel with meta information 
            id_channel = Channel.fromPath(file_integrity_glob).collect() // get all summary files
                .map{ it -> get_only_passing_samples(it)}.filter { it != null }.toList()
                 // loop through files and identify those that don't have "FAILED" in them and then parse file name and return those ids that pass

            // Get reads
            def r1_glob = append_to_path(params.indir.toString(),'*/fastp_trimd/*_1.trim.fastq.gz')
            def r2_glob = append_to_path(params.indir.toString(),'*/fastp_trimd/*_2.trim.fastq.gz')
            //create reads channel with meta information
            r1_reads_ch = Channel.fromPath(r1_glob).map{ it -> create_meta(it, "_1.trim.fastq.gz")}
            filtered_r1_reads_ch = r1_reads_ch.combine(id_channel).filter{ meta, read_1, id_channel -> id_channel.contains(meta.id)}
            .map{ meta, read_1, id_channel -> [ meta, read_1 ]}

            r2_reads_ch = Channel.fromPath(r2_glob).map{ it -> create_meta(it, "_2.trim.fastq.gz")}
            filtered_r2_reads_ch = r2_reads_ch.combine(id_channel).filter{ meta, read_2, id_channel -> id_channel.contains(meta.id)}
            .map{ meta, read_2, id_channel -> [ meta, read_2 ]}
            // combine reads into one channel
            combined_reads_ch = filtered_r1_reads_ch.join(filtered_r2_reads_ch, by: [0]).map{ meta, read_1, read_2 -> [meta, [read_1, read_2]]}

            // Get scaffolds
            def scaffolds_glob = append_to_path(params.indir.toString(),'*/assembly/*.filtered.scaffolds.fa.gz')
            //create scaffolds channel with meta information -- annoying, but you have to keep this in the brackets instead of having it once outside.
            scaffolds_ch = Channel.fromPath(scaffolds_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, ".filtered.scaffolds.fa.gz")} // create meta for sample
            // Checking regrex has correct extension
            scaffolds_ch = scaffolds_ch.map{ it -> check_scaffolds(it) }
            // Filter other channels based on meta.id
            filtered_scaffolds_ch = scaffolds_ch.combine(id_channel).filter{ meta, scaffolds, id_channel -> id_channel.contains(meta.id)}.map{ meta, scaffolds, id_channel -> [ meta, scaffolds ]}

            // get .tax files for MLST updating
            def taxa_glob = append_to_path(params.indir.toString(),'*/*.tax')
            //create .tax file channel with meta information 
            taxonomy_ch = Channel.fromPath(taxa_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, ".tax")} // create meta for sample
            //filtering out failured samples
            filtered_taxonomy_ch = taxonomy_ch.combine(id_channel).filter{ meta, taxa, id_channel -> id_channel.contains(meta.id)}.map{ meta, taxa, id_channel -> [ meta, taxa ]}

            // get prokka files
            def prokka_gff_glob = append_to_path(params.indir.toString(),'*/annotation/*.gff')
            def prokka_faa_glob = append_to_path(params.indir.toString(),'*/annotation/*.faa')
            //create .gff and .faa files channel with meta information 
            gff_ch = Channel.fromPath(prokka_gff_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, ".gff")} // create meta for sample
            //filtering out failured samples
            filtered_gff_ch = gff_ch.combine(id_channel).filter{ meta, gff, id_channel -> id_channel.contains(meta.id)}.map{ meta, gff, id_channel -> [ meta, gff ]}

            faa_ch = Channel.fromPath(prokka_faa_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, ".faa")} // create meta for sample
            //filtering out failured samples
            filtered_faa_ch = faa_ch.combine(id_channel).filter{ meta, faa, id_channel -> id_channel.contains(meta.id)}.map{ meta, faa, id_channel -> [ meta, faa ]}

            def line_summary_glob = append_to_path(params.indir.toString(),'*/*_summaryline.tsv')
            line_summary_ch = Channel.fromPath(line_summary_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, "_summaryline.tsv")} // create meta for sample
            //filtering out failured samples
            filtered_line_summary_ch = line_summary_ch.combine(id_channel).filter{ meta, line_summary, id_channel -> id_channel.contains(meta.id)}.map{ meta, line_summary, id_channel -> [ meta, line_summary ]}

            def synopsis_glob = append_to_path(params.indir.toString(),'*/*.synopsis')
            synopsis_ch = Channel.fromPath(synopsis_glob) // use created regrex to get samples
                .map{ it -> create_meta(it, ".synopsis")} // create meta for sample
            //filtering out failured samples
            filtered_synopsis_ch = synopsis_ch.combine(id_channel).filter{ meta, synopsis, id_channel -> id_channel.contains(meta.id)}.map{ meta, synopsis, id_channel -> [ meta, synopsis ]}

            def griphin_excel_glob = append_to_path(params.indir.toString(),'*_GRiPHin_Summary.xlsx')
            griphin_excel_ch = Channel.fromPath(griphin_excel_glob) // use created regrex to get samples
            def griphin_tsv_glob = append_to_path(params.indir.toString(),'*_GRiPHin_Summary.tsv')
            griphin_tsv_ch = Channel.fromPath(griphin_tsv_glob) // use created regrex to get samples

            //get valid samplesheet for griphin step in cdc_scaffolds
            /*CREATE_SAMPLESHEET (
                indir
            )
            ch_versions = ch_versions.mix(CREATE_SAMPLESHEET.out.versions)

            valid_samplesheet = CREATE_SAMPLESHEET.out.samplesheet*/

        } else if (samplesheet != null) {
            // if a samplesheet was passed then use that to create the channel

            SAMPLESHEET_CHECK (
                samplesheet, false, false, true
            )
            ch_versions = ch_versions.mix(SAMPLESHEET_CHECK.out.versions)

            directory_ch = SAMPLESHEET_CHECK.out.csv.splitCsv( header:true, sep:',' ).map{ create_dir_channels(it) }
            // get file_integrity file for MLST updating

            //create file_integrity file channel with meta information -- need to pass to DO_MLST subworkflow
            scaffolds_integrity_ch = directory_ch.map{dir -> append_to_path(dir,'file_integrity/', /.*_scaffolds_summary\.txt$/)} //use regrex to get files
            .map{ it -> create_meta(it, "_scaffolds_summary.txt")} // add meta information
            //scaffolds_integrity_ch.view()

            /*/get list of all samples in the folder - just using the file_integrity file, but it could be any file really
            file_integrity_ch = directory_ch.map{dir -> test_append_to_path(dir,'file_integrity/', /.*_\*_summary\.txt/)} //use regrex to get files
            file_integrity_ch.collect().view()
            //create file_integrity file channel with meta information 
            id_channel = file_integrity_ch.collect() // get all summary files
                .map{ it -> get_only_passing_samples(it)}.filter { it != null }.toList()
                // loop through files and identify those that don't have "FAILED" in them and then parse file name and return those ids that pass
            //id_channel.view()*/

            // Get reads
            r1_reads_ch = directory_ch.map{dir -> append_to_path(dir,'fastp_trimd/', /.*_1\.trim\.fastq\.gz$/)} //use regrex to get files
            .map{ it -> create_meta(it, "_1.trim.fastq.gz")} // add meta information
            r2_reads_ch = directory_ch.map{dir -> append_to_path(dir,'fastp_trimd/', /.*_2\.trim\.fastq\.gz$/)} //use regrex to get files
            .map{ it -> create_meta(it, "_2.trim.fastq.gz")} // add meta information
            /*create reads channel with meta information
            //filtered_r1_reads_ch = r1_reads_ch.combine(id_channel).filter{ meta, read_1, id_channel -> id_channel.contains(meta.id)}.map{ meta, read_1, id_channel -> [ meta, read_1 ]}
            // combine reads into one channel*/
            combined_reads_ch = r1_reads_ch.join(r2_reads_ch, by: [0]).map{ meta, read_1, read_2 -> [meta, [read_1, read_2]]}

            // Get scaffolds
            scaffolds_ch = directory_ch.map{dir -> append_to_path(dir,'assembly/', /.*\.filtered\.scaffolds\.fa\.gz$/)} //use regrex to get files
            .map{ it -> create_meta(it, ".filtered.scaffolds.fa.gz")} // add meta information
            .map{ it -> check_scaffolds(it) }// Checking regrex has correct extension
            //scaffolds_ch.view()
            //scaffolds_ch.map{ it -> check_scaffolds(it) }
            // Filter other channels based on meta.id
            //filtered_scaffolds_ch = scaffolds_ch.combine(id_channel).filter{ meta, scaffolds, id_channel -> id_channel.contains(meta.id)}.map{ meta, scaffolds, id_channel -> [ meta, scaffolds ]}
            //scaffolds_ch.view()

            /*/ get .tax files for MLST updating
            taxa_glob_ch = directory_ch.map{dir -> append_to_path(params.indir.toString(),'/*.tax')}
            //create .tax file channel with meta information 
            taxonomy_ch = taxa_glob_ch.flatMap{ regex -> Channel.fromPath(regex) } // use created regrex to get samples
                .map{ it -> create_meta(it, ".tax")} // create meta for sample
            //filtering out failured samples
            filtered_taxonomy_ch = taxonomy_ch.combine(id_channel).filter{ meta, taxa, id_channel -> id_channel.contains(meta.id)}.map{ meta, taxa, id_channel -> [ meta, taxa ]}

            // get prokka files
            prokka_gff_glob_ch = directory_ch.map{dir -> append_to_path(params.indir.toString(),'/annotation/*.gff$')}
            prokka_faa_glob_ch = directory_ch.map{dir -> append_to_path(params.indir.toString(),'/annotation/*.faa$')}
            //create .gff and .faa files channel with meta information 
            gff_ch = prokka_gff_glob_ch.flatMap{ regex -> Channel.fromPath(regex) } // use created regrex to get samples
                .map{ it -> create_meta(it, ".gff")} // create meta for sample
                //filtering out failured samples
            filtered_gff_ch = gff_ch.combine(id_channel).filter{ meta, gff, id_channel -> id_channel.contains(meta.id)}.map{ meta, gff, id_channel -> [ meta, gff ]}

            faa_ch = prokka_faa_glob_ch.flatMap{ regex -> Channel.fromPath(regex) } // use created regrex to get samples
                .map{ it -> create_meta(it, ".faa")} // create meta for sample
            //filtering out failured samples
            filtered_faa_ch = faa_ch.combine(id_channel).filter{ meta, faa, id_channel -> id_channel.contains(meta.id)}.map{ meta, faa, id_channel -> [ meta, faa ]}

            line_summary_glob_ch = directory_ch.map{dir -> append_to_path(params.indir.toString(),'/*_summaryline.tsv')}
            line_summary_ch = line_summary_glob_ch.flatMap{ regex -> Channel.fromPath(regex) } // use created regrex to get samples
                .map{ it -> create_meta(it, "_summaryline.tsv")} // create meta for sample
            //filtering out failured samples
            filtered_line_summary_ch = line_summary_ch.combine(id_channel).filter{ meta, line_summary, id_channel -> id_channel.contains(meta.id)}.map{ meta, line_summary, id_channel -> [ meta, line_summary ]}

            synopsis_glob_ch = directory_ch.map{dir -> append_to_path(params.indir.toString(),'/*.synopsis')}
            synopsis_ch = synopsis_glob_ch.flatMap{ regex -> Channel.fromPath(regex) } // use created regrex to get samples
                .map{ it -> create_meta(it, ".synopsis")} // create meta for sample
            //filtering out failured samples
            filtered_synopsis_ch = synopsis_ch.combine(id_channel).filter{ meta, synopsis, id_channel -> id_channel.contains(meta.id)}.map{ meta, synopsis, id_channel -> [ meta, synopsis ]}

            griphin_excel_glob_ch = directory_ch.map{dir -> append_to_path(params.indir.toString(),'*_GRiPHin_Summary.xlsx')}
            griphin_excel_ch = griphin_excel_glob_ch.flatMap{ regex -> Channel.fromPath(regex) } // use created regrex to get samples
            def griphin_tsv_glob = append_to_path(params.indir.toString(),'*_GRiPHin_Summary.tsv')
            griphin_tsv_ch = Channel.fromPath(griphin_tsv_glob) // use created regrex to get samples*/

            valid_samplesheet = SAMPLESHEET_CHECK.out.csv

        } else {
            exit 1, 'You need EITHER an input directory samplesheet (using --input) or a single sample directory (using --indir)!' 
        }

    emit:
        filtered_scaffolds = scaffolds_ch      // channel: [ meta, [ scaffolds_file ] ]
        reads              = combined_reads_ch
       /* taxonomy           = filtered_taxonomy_ch
        prokka_gff         = filtered_gff_ch
        prokka_faa         = filtered_faa_ch*/
        fairy_outcome      = scaffolds_integrity_ch
        /*line_summary       = filtered_line_summary_ch
        synopsis           = filtered_synopsis_ch*/
        valid_samplesheet  = valid_samplesheet
        versions           = ch_versions

}

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/


// Function to get list of [ meta, [ directory ] ]
def create_dir_channels(LinkedHashMap row) {
    array = [ row.directory ]
    return array
}

/*/ Function to filter files and parse their names
def get_only_passing_samples(filePath) {
    filePaths.each { filePath ->
        def fileObj = new File(filePath)
        if (!fileObj.text.contains('FAILED')) {
            // Get the file name by removing the path and extension
            def fileName = filePath.getFileName().toString()
            println(fileName)
            // Remove everything after the regex '_*_summary.txt'
            return fileName.replaceAll(/_.+_summary\.txt/, '')
        }
        return null // Return null for files containing 'FAILED'
    }
}*/

// Define a function to append glob patterns and create meta
def append_and_create_meta(dir, type) {
    if (type = "file_integrity") {
        def globPattern = "${dir}/*/file_integrity/*_scaffolds_summary.txt"
        return Channel.fromPath(globPattern).map { path -> create_meta(path, "_scaffolds_summary.txt")
        }
    }
}

/*/ Function to filter files and parse their names
def get_only_passing_samples(files) {
    def passingFiles = []
    files.each { file ->
        def fileObj = new File(file.toString())
        if (!fileObj.text.contains('FAILED')) {
            // Get the file name by removing everything before the last '/'
            def fileName = file.toString().tokenize('/').last()
            // Remove everything after the regex '_*_summary.txt'
            def parsedName = fileName.replaceAll(/_.*_summary\.txt/, '')
            //println "Parsed file name: ${parsedName}"
            passingFiles << parsedName
        }
    }
    return passingFiles
}*/

def append_to_path(partial_path, phx_path, regrex) {
    if (partial_path[0].toString().endsWith('/')) {
        //new_string = partial_path[0].toString() + phx_path + regrex
        ab_path = partial_path[0].toString() + phx_path
    } else {
        //new_string = partial_path[0].toString()+ '/' + phx_path + regrex
        ab_path = partial_path[0].toString() + '/' + phx_path
    }
    // get files 
    def fileList = []
    def dir = new File(ab_path)
    if (dir.exists() && dir.isDirectory()) {
        dir.eachFileRecurse { file ->
            if (file.isFile() && file.name ==~ regrex) {
                fileList << file.absolutePath
            }
        }
    }
    return fileList
}

def test_append_to_path(partial_path, phx_path, regrex) {
    if (partial_path[0].toString().endsWith('/')) {
        new_string = partial_path[0].toString() + phx_path + regrex
        ab_path = partial_path[0].toString() + phx_path
    } else {
        new_string = partial_path[0].toString()+ '/' + phx_path + regrex
        ab_path = partial_path[0].toString() + '/' + phx_path
    }
    // get files 
    def fileList = []
    def dir = new File(ab_path)
    if (dir.exists() && dir.isDirectory()) {
        dir.eachFileRecurse { file ->
            if (file.isFile() && file.name ==~ regrex) {
                fileList << file.absolutePath
            }
        }
    }
    return fileList
}

def create_meta(sample, file_extension){
    '''Creating meta: [[id:sample1], $PATH/sample1.filtered.scaffolds.fa.gz]'''
    sample_name_minus_path = sample[0].toString().split('/')[-1] // get the last string after the last backslash
    sample_name = sample_name_minus_path.replaceAll(file_extension, "") // remove file extention to get only sample name
    def meta = [:] // create meta array
    meta.id = sample_name
    array = [ meta, sample ]  //file() portion provides full path
    return array
}

def check_scaffolds(scaffold_channel) {
    if (scaffold_channel[1][0].toString().endsWith(".fasta.gz") | scaffold_channel[1][0].toString().endsWith(".fa.gz") ) {
        //println(scaffold_channel)
        //println(scaffold_channel[1][0])
        return scaffold_channel
    } else {
        println(scaffold_channel)
        println(scaffold_channel[1][0])
        exit 1, "ERROR: No scaffolds found. Either your scaffold regrex is off (scaffolds files should end in either '.fa.gz' or ''.fasta.gz') or the directory provided doesn't contain assembly files." 
    }
}