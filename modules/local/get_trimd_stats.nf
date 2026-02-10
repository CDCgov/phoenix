process GET_TRIMD_STATS {
    tag "$meta.id"
    label 'process_single'
    stageInMode 'copy'
    // base_v2.2.0 - MUST manually change below (line 30)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta), path(fastp_trimd_json), path(fastp_singles_json), path(raw_qc), path(fairy_outcome)
    val(busco_val)
    val(phx_version)

    output:
    tuple val(meta), path('*_trimmed_read_counts.txt'),               emit: fastp_total_qc
    tuple val(meta), path('*_trimstats_summary.txt'),  optional:true, emit: outcome
    path('*_summaryline.tsv'),                         optional:true, emit: summary_line
    tuple val(meta), path('*_summary_old_3.txt'),                     emit: outcome_to_edit
    tuple val(meta), path('*.synopsis'),               optional:true, emit: synopsis
    path("versions.yml"),                                             emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def busco_parameter = busco_val ? "--busco" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}FastP_QC.py \\
        --trimmed_json ${fastp_trimd_json} \\
        --single_json ${fastp_singles_json} \\
        --name ${prefix}

    # Check that there are still reads in R1 and R2 before fastqc. If there aren't reads then fastqc dies.

    # Output check for messages indicating there are no trimmed reads after filtering.
    ${ica}fairy.py -r ${raw_qc} -f ${fairy_outcome} -t ${prefix}_trimmed_read_counts.txt ${busco_parameter} --phx_version $phx_version

    #making a copy of the summary file to pass to BBMAP_REFORMAT to handle file names being the same
    cp ${prefix}_trimstats_summary.txt ${prefix}_summary_old_3.txt
    if ! grep -q "FAILED" ${prefix}_trimstats_summary.txt; then
        rm ${prefix}_trimstats_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        fairy.py: \$( ${ica}fairy.py --version )
        FastP_QC.py: \$(${ica}FastP_QC.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}