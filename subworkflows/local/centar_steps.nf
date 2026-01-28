/*
========================================================================================
    IMPORT LOCAL MODULES
========================================================================================
*/

include { CDIFF_CLADE                    } from '../../modules/local/centar/cdiff_clade'
include { CDIFF_TOXINOTYPER              } from '../../modules/local/centar/cdiff_toxinotyper'
//include { CDIFF_PLASMIDS                 } from '../../modules/local/centar/cdiff_plasmids'
include { GAMMA as CDIFF_TOX_GENES       } from '../../modules/local/gamma'
//include { GAMMA as CDIFF_AR_GENES_ALL    } from '../../modules/local/gamma'
include { GAMMA as CDIFF_AR_GENES_AA     } from '../../modules/local/gamma'
include { GAMMA_S as CDIFF_AR_GENES_NT   } from '../../modules/local/gammas'
include { WGMLST                         } from '../../modules/local/centar/wgmlst'
include { CDIFF_RIBOTYPER                } from '../../modules/local/centar/cdiff_ribotyper'
include { CENTAR_CONSOLIDATER            } from '../../modules/local/centar/centar_consolidater'
include { KRAKEN2_WF as KRAKEN2_PLASMID  } from './kraken2krona'

/*
========================================================================================
    GROOVY FUNCTIONS
========================================================================================
*/

// Groovy funtion to make [ meta.id, [] ] - just an empty channel
def create_empty_ch(input_for_meta) { // We need meta.id associated with the empty list which is why .ifempty([]) won't work
    meta_id = input_for_meta[0]
    output_array = [ meta_id, [] ]
    return output_array
}

def get_taxa(input_ch){ 
        def genus = ""
        def species = ""
        input_ch[1].eachLine { line ->
            try {
                if (line.startsWith("G:")) {
                    // Check if we have the format with ID number like "G:32008 Burkholderia"
                    String remaining = line.substring(2).trim();
                    
                    // Check if the remaining part starts with a number followed by space
                    if (remaining.matches("^\\d+\\s+.*")) {
                        // For format like "G:32008 Burkholderia"
                        genus = remaining.replaceFirst("^\\d+\\s+", "").trim();
                    } else {
                        // For format like "G:      Serratia"
                        genus = remaining;
                    }
                } else if (line.startsWith("s:")) {
                    // Same logic for species
                    String remaining = line.substring(2).trim();
                    
                    if (remaining.matches("^\\d+\\s+.*")) {
                        // For format like "s:95486 cenocepacia"
                        species = remaining.replaceFirst("^\\d+\\s+", "").trim();
                    } else {
                        // For format like "s:      marcescens"
                        species = remaining;
                    }
                }
            } catch (ArrayIndexOutOfBoundsException e) {
                // Handle specific "Index 1 out of bounds for length 1" error - when there are is only a genus. line.split(":")[1].trim().split('\t')[1] - [1] is not needed
                if (e.getMessage() != null) {
                    if (e.message.contains("Index 1 out of bounds for length 1")) {
                        if (line.startsWith("G:")) {
                            // Use the fallback logic if accessing index 1 fails
                            genus = line.split(":")[1].trim()
                        } else if (line.startsWith("s:")) {
                            species = line.split(":")[1].trim()
                        }
                    } else {
                        println("E:" + e + "see get_taxa() in centar_steps.nf, parsing of *.tax failed.")
                        // Re-throw or handle other ArrayIndexOutOfBoundsExceptions if necessary
                        println "Unexpected ArrayIndexOutOfBoundsException: ${e.message}"
                    }
                } else {
                    println("Error message is empty?!")
                }

                
            } catch (IndexOutOfBoundsException e) {
                // Handle specific "Index 1 out of bounds for length 1" error - when there are is only a genus. line.split(":")[1].trim().split('\t')[1] - [1] is not needed
                if (e.getMessage() != null) {
                    if (e.message.contains("Index 1 out of bounds for length 1")) {
                        if (line.startsWith("G:")) {
                            // Use the fallback logic if accessing index 1 fails
                            genus = line.split(":")[1].trim()
                        } else if (line.startsWith("s:")) {
                            species = line.split(":")[1].trim()
                        }
                    } else {
                        println("E:" + e + "see get_taxa() in centar_steps.nf, parsing of *.tax failed.")
                        // Re-throw or handle other IndexOutOfBoundsExceptions if necessary
                        println "Unexpected IndexOutOfBoundsException: ${e.message}"
                    }
                } else {
                    println("Error message is empty?!")
                }
            } catch (Exception e) {
                // Catch any other exceptions that might occur
                println "An unexpected error occurred in get_taxa() in centar_steps.nf: ${e.message}"
            }
        }
    // Change from Genus species match, to just species. Since Clostridioidesd only has a very small list of species we will go ahead and check all to gauge toxicity
    //return [input_ch[0], "$genus $species" ]
    return [ input_ch[0], "$genus" ]
}

def get_only_taxa(input_ch){ 
        def genus = ""
        def species = ""
        input_ch[1].eachLine { line ->
            if (line.startsWith("G:")) {
                genus = line.split(":")[1].trim().split('\t')[1]
            } else if (line.startsWith("s:")) {
                species = line.split(":")[1].trim().split('\t')[1]
            }
        }
        return [ "$genus $species" ]
}

def check_params_var(species_bol, species_param) {
    if (species_bol == true && species_param == true){
        return true
    } else if (species_bol == true && species_param == false) {
        return false
    } else {
        return false
    }
}

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow CENTAR_SUBWORKFLOW {
    take:
        combined_mlst       // CREATE_INPUT_CHANNELS.out.combined_mlst
        fairy_outcome       // CREATE_INPUT_CHANNELS.out.fairy_outcome
        filtered_scaffolds  // CREATE_INPUT_CHANNELS.out.filtered_scaffolds
        mlst_db             // ASSET_CHECK.out.mlst_db
        taxonomy            // taxa information

    main:
        // Allow outdir to be relative
        outdir_path = Channel.fromPath(params.outdir, relative: true)
        ch_versions = Channel.empty() // Used to collect the software versions

        // Check the taxa file to confirm the organism is Clostridioides difficile
        clade_ch = combined_mlst.join(taxonomy.map{it -> get_taxa(it)}, by: [[0][0],[0][1]])\
                    .filter{meta, combined_mlst, cdiff_check -> cdiff_check.tokenize()[0] == "Clostridioides" }\
                    .map{   meta, combined_mlst, cdiff_check -> [meta, combined_mlst] }.combine(mlst_db)

        // get c diff clade
        CDIFF_CLADE (
            clade_ch
        )
        ch_versions = ch_versions.mix(CDIFF_CLADE.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        cdiff_check = taxonomy.map{it -> get_taxa(it)} // get organism from file
        // we will filter out an samples that aren't Clostridioides difficile
        filtered_scaffolds_ch = filtered_scaffolds.map{        meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}\
            .join(fairy_outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome      -> [[id:meta.id, project_id:meta.project_id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])\
            .filter {                                          meta, filtered_scaffolds, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}\
            .map{                                              meta, filtered_scaffolds, fairy_outcome -> return [meta, filtered_scaffolds] }\
            .join(cdiff_check, by: [[0][0],[0][1]]).filter{    meta, filtered_scaffolds, taxa_id -> taxa_id.tokenize()[0] == "Clostridioides" }.map{ meta, filtered_scaffolds, taxa_id -> [meta, filtered_scaffolds ]}

        // Running gamma to identify toxin genes in scaffolds for general presence
        CDIFF_TOX_GENES (
            filtered_scaffolds_ch, params.cdiff_tox_gene_db
        )
        ch_versions = ch_versions.mix(CDIFF_TOX_GENES.out.versions)

        /*/ Running gamma to identify Cdiff specific AR genes in scaffolds
        CDIFF_AR_GENES_ALL (
            filtered_scaffolds_ch, params.cdiff_ar_gene_ALL_db
        )
        ch_versions = ch_versions.mix(CDIFF_AR_GENES_ALL.out.versions)*/

        // Running gamma to identify Cdiff specific AR genes in scaffolds
        CDIFF_AR_GENES_AA (
            filtered_scaffolds_ch, params.cdiff_ar_gene_AA_db
        )
        ch_versions = ch_versions.mix(CDIFF_AR_GENES_AA.out.versions)

        // Running gamma to identify Cdiff specific AR genes in scaffolds
        CDIFF_AR_GENES_NT (
            filtered_scaffolds_ch, params.cdiff_ar_gene_NT_db
        )
        ch_versions = ch_versions.mix(CDIFF_AR_GENES_NT.out.versions)

        // Running blat to identify diffbase toxin genes for specific toxinotyping
        CDIFF_TOXINOTYPER (
            filtered_scaffolds_ch, params.cdiff_diffbase_AA, params.cdiff_diffbase_definitions
        )
        ch_versions = ch_versions.mix(CDIFF_TOXINOTYPER.out.versions)

        if (params.wgmlst_container != null) {
            // Running blat to identify diffbase toxin genes for specific toxinotyping
            WGMLST (
                filtered_scaffolds_ch, params.blast_kb, params.cdiff_wgmlst_blast_db, params.cdiff_wgmlst_blast_loci, params.qckb
            )
            ch_versions = ch_versions.mix(WGMLST.out.versions)

            // Join everything together based on meta.id and project_id
            allele_calls_ch = WGMLST.out.csv_core.map{meta, csv_core_st, csv_core_pcr           -> [[id:meta.id, project_id:meta.project_id], csv_core_st]}\
                .join(WGMLST.out.csv_accessory.map{   meta, csv_accessory_st, csv_accessory_pcr -> [[id:meta.id, project_id:meta.project_id], csv_accessory_st]}, by: [[0][0],[0][1]])

            // Running blat to identify diffbase toxin genes for specific toxinotyping
            CDIFF_RIBOTYPER (
                allele_calls_ch.map{meta, csv_core, csv_accessory -> [ meta, csv_core ]},
                allele_calls_ch.map{meta, csv_core, csv_accessory -> [ meta, csv_accessory ]}
            )
            ch_versions = ch_versions.mix(CDIFF_RIBOTYPER.out.versions)

            detailed_ribotype_file_ch = CDIFF_RIBOTYPER.out.detailed_ribotype_file

        } else {
            //create an empty channel for ribotyper if it is not run -- this causes a failure right now
            detailed_ribotype_file_ch = combined_mlst.map{ it -> create_empty_ch(it)}
        }

        // Join everything together based on meta.id
        cdiff_summary_ch = CDIFF_TOX_GENES.out.gamma.map{meta, gamma                    -> [[id:meta.id, project_id:meta.project_id], gamma]}\
            .join(CDIFF_CLADE.out.clade.map{             meta, clade                    -> [[id:meta.id, project_id:meta.project_id], clade]},                  by: [[0][0],[0][1]])\
            .join(CDIFF_TOXINOTYPER.out.tox_file.map{    meta, tox_file                 -> [[id:meta.id, project_id:meta.project_id], tox_file]},               by: [[0][0],[0][1]])\
            .join(CDIFF_AR_GENES_AA.out.gamma.map{       meta, gamma                    -> [[id:meta.id, project_id:meta.project_id], gamma]},                  by: [[0][0],[0][1]])\
            .join(CDIFF_AR_GENES_NT.out.gamma.map{       meta, gamma                    -> [[id:meta.id, project_id:meta.project_id], gamma]},                  by: [[0][0],[0][1]])\
            .join(detailed_ribotype_file_ch.map{         meta, detailed_ribotype_file   -> [[id:meta.id, project_id:meta.project_id], detailed_ribotype_file]}, by: [[0][0],[0][1]])

        CENTAR_CONSOLIDATER (
            cdiff_summary_ch, params.cemb_strt_xwalk
        )
        ch_versions = ch_versions.mix(CENTAR_CONSOLIDATER.out.versions)

    emit:
        consolidated_centar = CENTAR_CONSOLIDATER.out.centar_summary_line
        versions            = ch_versions // channel: [ versions.yml ]

}
