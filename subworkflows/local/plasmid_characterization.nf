//
// Subworkflow: run blast and ani and VF gamma on plasmid contigs 
//


include { PLSDB       } from '../../modules/local/long_read/plsdb'   //module to blast contigs against plsdb
include { PLSDBANI  } from '../../modules/local/long_read/plsdbani' //module to run fastani on contigs against plsdb accessions
include { GAMMA_LR  } from '../../modules/local/long_read/gamma_LR' //module to run gamma with VF db

workflow PLASMID_CHARACTERIZATION {
    take:    //declare input channels
        fasta     // channel: tuple val(meta), path(fasta)
        //scaffold_count_check // SCAFFOLD_COUNT_CHECK.out.outcome
        plsdb      // channel: path(plsdb_dir)
        conf   // channel: path(conf)
        viz     // channel: path(viz)
        plsdbfasta // channel: path(plsdbfasta)
        vfdb    // channel: path(vfdb) - virulence factor database for plasmids
    main:
        ch_versions = Channel.empty() // Used to collect the software versions
        // Creating channel to ensure ID is paired with matching fasta
        fasta_ch = fasta.map{meta, fasta -> [[id:meta.id], fasta]}\
        //.join(scaffold_count_check.splitCsv(strip:true, by:5).map{meta, fairy_outcome -> [[id:meta.id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [0])\
        //.filter { meta, reads, fairy_outcome -> fairy_outcome[4] == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."}\
        //.map{ meta, fasta, fairy_outcome -> return [meta, fasta] }

        // Running blast against plsdb
        PLSDB (      
            fasta_ch, plsdb, conf       
        )
        ch_versions = ch_versions.mix(PLSDB.out.versions)  
        
        // Running fastani to compare plsdb hits to contigs
        plsdb_ch = PLSDB.out.plasmidID.map{meta, plasmidID -> [[id:meta.id], plasmidID]}\
        .join(fasta_ch, by: [0])\
        .map{meta, plasmidID, fasta -> [[id:meta.id], plasmidID, fasta]}
                
        PLSDBANI (
            plsdb_ch, viz, plsdbfasta
        )
        ch_versions = ch_versions.mix(PLSDBANI.out.versions)   

        // Running gamma to find virulence factors in plasmid contigs
        GAMMA_LR (
            fasta_ch, vfdb
        )
        ch_versions = ch_versions.mix(GAMMA_LR.out.versions)

    emit:
        blast_out = PLSDB.out.plasmidID
        ani_out   = PLSDBANI.out.plasmidANI
        plasmid_vf_out = GAMMA_LR.out.vf
        versions  = ch_versions // channel: [ versions.yml ]        


}