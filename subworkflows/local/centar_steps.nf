/*
========================================================================================
    IMPORT LOCAL MODULES
========================================================================================
*/

include { CDIFF_CLADE                    } from '../../modules/local/centar/cdiff_clade'
include { CDIFF_TOXINOTYPER              } from '../../modules/local/centar/cdiff_toxinotyper'
include { CDIFF_PLASMIDS                 } from '../../modules/local/centar/cdiff_plasmids'
include { GAMMA as CDIFF_TOX_GENES       } from '../../modules/local/gamma'
include { GAMMA as CDIFF_AR_GENES        } from '../../modules/local/gamma'
include { WGMLST                         } from '../../modules/local/centar/wgmlst'
include { CDIFF_RIBOTYPER                } from '../../modules/local/centar/cdiff_ribotyper'
include { CENTAR_CONSOLIDATER            } from '../../modules/local/centar/centar_consolidater'
include { KRAKEN2_WF as KRAKEN2_PLASMID  } from './kraken2krona'

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

    main:
        // Allow outdir to be relative
        outdir_path = Channel.fromPath(params.outdir, relative: true)
        ch_versions = Channel.empty() // Used to collect the software versions

        CDIFF_CLADE (
            combined_mlst, mlst_db
        )
        ch_versions = ch_versions.mix(CDIFF_CLADE.out.versions)

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        filtered_scaffolds_ch = filtered_scaffolds.map{    meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}
        .join(fairy_outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome      -> [[id:meta.id, project_id:meta.project_id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])

        CDIFF_PLASMIDS (
            filtered_scaffolds_ch, params.cdiff_plasmid_db
        )
        ch_versions = ch_versions.mix(CDIFF_PLASMIDS.out.versions)

        // Running gamma to identify toxin genes in scaffolds for general presence
        CDIFF_TOX_GENES (
            filtered_scaffolds_ch, params.cdiff_tox_gene_db
        )
        ch_versions = ch_versions.mix(CDIFF_TOX_GENES.out.versions)

        // Running gamma to identify Cdiff specific AR genes in scaffolds
        CDIFF_AR_GENES (
            filtered_scaffolds_ch, params.cdiff_ar_gene_db
        )
        ch_versions = ch_versions.mix(CDIFF_AR_GENES.out.versions)

        // Running blat to identify diffbase toxin genes for specific toxinotyping
        CDIFF_TOXINOTYPER (
            filtered_scaffolds_ch, params.cdiff_diffbase_AA, params.cdiff_diffbase_definitions
        )
        ch_versions = ch_versions.mix(CDIFF_TOXINOTYPER.out.versions)

        // Running blat to identify diffbase toxin genes for specific toxinotyping
        WGMLST (
            filtered_scaffolds_ch, params.cdiff_wgmlst_blast_db
        )
        ch_versions = ch_versions.mix(WGMLST.out.versions)

        // Running blat to identify diffbase toxin genes for specific toxinotyping
        CDIFF_RIBOTYPER (
            WGMLST.out.wgmlst_alleles_file
        )
        ch_versions = ch_versions.mix(CDIFF_RIBOTYPER.out.versions)

        // Join everything together based on meta.id
        cdiff_summary_ch = CDIFF_TOX_GENES.out.gamma.map{        meta, gamma           -> [[id:meta.id], gamma]}\
        .join(CDIFF_CLADE.out.clade.map{                         meta, clade           -> [[id:meta.id], clade]},         by: [0])\
        .join(CDIFF_TOXINOTYPER.out.tox_file.map{                meta, tox_file        -> [[id:meta.id], tox_file]},      by: [0])\
        .join(CDIFF_AR_GENES.out.gamma.map{                      meta, gamma           -> [[id:meta.id], gamma]},         by: [0])\
        .join(CDIFF_PLASMIDS.out.plasmids_file.map{              meta, plasmids_file   -> [[id:meta.id], plasmids_file]}, by: [0])\
        .join(CDIFF_RIBOTYPER.out.ribotype_file.map{             meta, ribotype_file   -> [[id:meta.id], ribotype_file]}, by: [0])

        CENTAR_CONSOLIDATER(
           cdiff_summary_ch
        )
        ch_versions = ch_versions.mix(CENTAR_CONSOLIDATER.out.versions)

    emit:
        consolidated_centar = CENTAR_CONSOLIDATER.out.centar_summary_line
        versions            = ch_versions // channel: [ versions.yml ]

}
