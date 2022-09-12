process CALCULATE_ASSEMBLY_RATIO {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.0.0'

    input:
    tuple val(meta), path(taxa_file), path(quast_report)
    path(ncbi_database)

    output:
    tuple val(meta), path('*_Assembly_ratio_*.txt'), emit: ratio
    tuple val(meta), path('*_GC_content_*.txt')     , emit: gc_content

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    calculate_assembly_ratio.sh -d $ncbi_database -q $quast_report -t $taxa_file -s ${prefix}
    """
}
