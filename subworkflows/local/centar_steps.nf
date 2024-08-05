/*
========================================================================================
    IMPORT LOCAL MODULES
========================================================================================
*/

include { CDIFF_CLADE                    } from '../../modules/local/centar/cdiff_clade'
include { CDIFF_PLASMID                  } from '../../modules/local/centar/cdiff_plasmid'
include { GAMMA as CDIFF_GENES           } from '../../modules/local/gamma'
include { CENTAR_CONSOLIDATER            } from '../../modules/local/centar/centar_consolidater'

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

        CDIFF_CLADE(
            combined_mlst, mlst_db
        )
        ch_versions = ch_versions.mix(CDIFF_CLADE.out.versions)

        /*CDIFF_PLASMID( // not written yet, just a placeholder
            ????
        )
        ch_versions = ch_versions.mix(CDIFF_PLASMID.out.versions)*/

        //combing scaffolds with scaffold check information to ensure processes that need scaffolds only run when there are scaffolds in the file
        filtered_scaffolds_ch = filtered_scaffolds.map{    meta, filtered_scaffolds -> [[id:meta.id, project_id:meta.project_id], filtered_scaffolds]}
        .join(fairy_outcome.splitCsv(strip:true, by:5).map{meta, fairy_outcome      -> [[id:meta.id, project_id:meta.project_id], [fairy_outcome[0][0], fairy_outcome[1][0], fairy_outcome[2][0], fairy_outcome[3][0], fairy_outcome[4][0]]]}, by: [[0][0],[0][1]])

        // Running gamma to identify hypervirulence genes in scaffolds
        CDIFF_GENES (
            filtered_scaffolds_ch, params.cdiff_gene_db
        )
        ch_versions = ch_versions.mix(CDIFF_GENES.out.versions)

        CENTAR_CONSOLIDATER(
           CDIFF_GENES.out.gamma, CDIFF_CLADE.out.clade
        )
        ch_versions = ch_versions.mix(CENTAR_CONSOLIDATER.out.versions)

    emit:
        consolidated_centar = CENTAR_CONSOLIDATER.out.report
        versions            = ch_versions // channel: [ versions.yml ]

}