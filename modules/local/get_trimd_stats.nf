process GET_TRIMD_STATS {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 30)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(fastp_trimd_json),
    path(fastp_singles_json),
    path(raw_qc), 
    path(fairy_outcome)
    val(busco_val)

    output:
    tuple val(meta), path('*_trimmed_read_counts.txt'),          emit: fastp_total_qc
    path('*_summaryline.tsv'),                                   optional:true, emit: summary_line
    tuple val(meta), path('*_summary_fastp.txt'),                emit: outcome
    tuple val(meta), path('*.synopsis'),                         optional:true, emit: synopsis
    path("versions.yml"),                                        emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_parameter = busco_val ? "--busco" : ""
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script_fastp = params.ica ? "python ${params.ica_path}/FastP_QC.py" : "FastP_QC.py"
    def script_fairy = params.ica ? "python ${params.ica_path}/fairy.py" : "fairy.py"
    """
    ${script_fastp} \\
      --trimmed_json ${fastp_trimd_json} \\
      --single_json ${fastp_singles_json} \\
      --name ${prefix}

    # making a copy of the summary file - this avoids writing to the previous file
    cp ${fairy_outcome} ${prefix}_input.txt

    # Output check for messages indicating there are no trimmed reads after filtering.
    ${script_fairy} -r ${raw_qc} -f ${prefix}_input.txt -t ${prefix}_trimmed_read_counts.txt ${busco_parameter}

    #making a copy of the summary file to pass to BBMAP_REFORMAT to handle file names being the same
    mv ${prefix}_summary.txt ${prefix}_summary_fastp.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        fairy.py: \$( ${script_fairy} --version )
        FastP_QC.py: \$(${script_fastp} --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}