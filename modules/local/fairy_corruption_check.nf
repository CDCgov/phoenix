process CORRUPTION_CHECK {
    tag "${meta.id}"
    label 'process_medium'
    // base_v2.1.0 - MUST manually change below (line 28)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(reads)
    val(busco_val)

    output:
    tuple val(meta), path('*_summary.txt'),                    emit: outcome
    tuple val(meta), path('*_summary_old.txt'),                emit: outcome_to_edit
    path('*_summaryline.tsv'),                  optional:true, emit: summary_line
    tuple val(meta), path('*.synopsis'),        optional:true, emit: synopsis
    path("versions.yml"),                                      emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def num1 = "${reads[0]}".minus(".fastq.gz")
    def num2 = "${reads[1]}".minus(".fastq.gz")
    def busco_parameter = busco_val ? "-b" : ""
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    #set +e
    #check for file integrity and log errors
    #if there is a corruption problem the script will create a *_summaryline.tsv and *.synopsis file for the sample.
    ${ica}fairy_proc.sh -r ${reads[0]} -p ${prefix} ${busco_parameter}
    ${ica}fairy_proc.sh -r ${reads[1]} -p ${prefix} ${busco_parameter}

    script_version=\$(${ica}fairy_proc.sh -V)

    #making a copy of the summary file to pass to READ_COUNT_CHECKS to handle file names being the same
    cp ${prefix}_summary.txt ${prefix}_summary_old.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}