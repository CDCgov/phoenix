process CALCULATE_ASSEMBLY_RATIO {
    tag "$meta.id"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 27)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta), path(taxa_file), path(quast_report)
    path(ncbi_database)

    output:
    tuple val(meta), path('*_Assembly_ratio_*.txt'), emit: ratio
    tuple val(meta), path('*_GC_content_*.txt')    , emit: gc_content
    path("versions.yml")                           , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // terra=true sets paths for bc/wget for terra container paths
    def terra = params.terra ? "-t terra" : ""
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "bash ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}calculate_assembly_ratio.sh -d $ncbi_database -q $quast_report -x $taxa_file -s ${prefix} $terra

    script_version=\$(${ica}calculate_assembly_ratio.sh -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NCBI_Assembly_Stats_DB: $ncbi_database
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}
