process GET_RAW_STATS {
    tag "${meta.id}"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(reads), val(fairy_corrupt_outcome)

    output:
    tuple val(meta), path('*_stats.txt'),           emit: raw_stats
    tuple val(meta), path('*_raw_read_counts.txt'), emit: combined_raw_stats
    path("versions.yml"),                           emit: versions

    when:
    //if the files are not corrupt then get the read stats
    "${fairy_corrupt_outcome[0]}" == "PASSED: File ${meta.id}_R1 is not corrupt." && "${fairy_corrupt_outcome[1]}" == "PASSED: File ${meta.id}_R2 is not corrupt."

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}q30.py ${reads[0]} > ${prefix}_R1_stats.txt
    ${ica}q30.py ${reads[1]} > ${prefix}_R2_stats.txt
    ${ica}create_raw_stats_output.py -n ${prefix} -r1 ${prefix}_R1_stats.txt -r2 ${prefix}_R2_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}