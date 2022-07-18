process DETERMINE_TAXA_ID_FAILURE {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(k2_bh_summary), val(spades_outcome)
    path(taxa_file)

    output:
    tuple val(meta), path('*.tax'), emit: taxonomy

    when:
    "${spades_outcome[0]}" == "run_failure" || "${spades_outcome[1]}" == "no_scaffolds" || "${spades_outcome[2]}" == "no_contigs"

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    determine_taxID.sh -r $k2_bh_summary -s $meta.id -d $taxa_file 
    """
}