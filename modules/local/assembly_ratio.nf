process CALCULATE_ASSEMBLY_RATIO {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.0.2'

    input:
    tuple val(meta), path(taxa_file), path(quast_report)
    path(ncbi_database)

    output:
    tuple val(meta), path('*_Assembly_ratio_*.txt'), emit: ratio
    tuple val(meta), path('*_GC_content_*.txt')    , emit: gc_content
    path("versions.yml")                           , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) {
        terra = ""
    } else if (params.terra==true) {
        terra = "-t terra"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    calculate_assembly_ratio.sh -d $ncbi_database -q $quast_report -x $taxa_file -s ${prefix} $terra

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NCBI Assembly Stats DB: $ncbi_database
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
