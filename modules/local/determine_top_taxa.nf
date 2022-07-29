process DETERMINE_TOP_TAXA {
    tag "$meta.id"
    label 'process_low'
    //container 'staphb/gamma:2.1'

    input:
    tuple val(meta), path(mash_dists), path(assembly_scaffolds)
    path(refseq_fasta_database)

    output:
    tuple val(meta), path('*_best_MASH_hits.txt')     , emit: top_taxa_list
    tuple val(meta), path('Ref_Seq/*_genomic.fna.gz') , emit: reference_files

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sort_and_prep_dist.sh -a $assembly_scaffolds -x $mash_dists -d $refseq_fasta_database

    rsync -r Ref_Seq "${baseDir}/assets/databases/"
    """
}