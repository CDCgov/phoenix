process GET_RAW_STATS {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 32)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(reads), path(fairy_outcome), val(fairy_corrupt_outcome)
    val(busco_val)

    output:
    tuple val(meta), path('*_stats.txt'),                        emit: raw_stats
    tuple val(meta), path('*_raw_read_counts.txt'),              emit: combined_raw_stats
    tuple val(meta), path('*_summary_rawstats.txt'),             emit: outcome
    path('*_summaryline.tsv'),                                   optional:true, emit: summary_line
    tuple val(meta), path('*.synopsis'),                         optional:true, emit: synopsis
    path("versions.yml"),                                        emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_parameter = busco_val ? "--busco" : ""
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script_q30 = params.ica ? "python ${params.ica_path}/q30.py" : "q30.py"
    def script_stats = params.ica ? "python ${params.ica_path}/create_raw_stats_output.py" : "create_raw_stats_output.py"
    def script_fairy = params.ica ? "python ${params.ica_path}/fairy.py" : "fairy.py"
    """
    ${script_q30} -i ${reads[0]} > ${prefix}_R1_stats.txt
    ${script_q30} -i ${reads[1]} > ${prefix}_R2_stats.txt
    ${script_stats} -n ${prefix} -r1 ${prefix}_R1_stats.txt -r2 ${prefix}_R2_stats.txt

    # making a copy of the summary file - this avoids writing to the previous file
    cp ${fairy_outcome} ${prefix}_input.txt

    # Output check for messages indicating read pairs that do not match
    ${script_fairy} -r ${prefix}_raw_read_counts.txt -f ${prefix}_input.txt ${busco_parameter}

    # rename output file
    mv ${prefix}_summary.txt ${prefix}_summary_rawstats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        q30.py: \$(${script_q30} --version )
        create_raw_stats_output.py: \$(${script_stats} --version )
        fairy.py: \$(${script_fairy} --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}