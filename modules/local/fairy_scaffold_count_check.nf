process SCAFFOLD_COUNT_CHECK {
    tag "${meta.id}"
    label 'process_medium'
    // base_v2.1.0 - MUST manually change below (line 50)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(bbmap_log), path(fairy_read_count_outcome),
    path(raw_qc_file),
    path(fastp_total_qc_file),
    path(kraken2_trimd_report_file),
    path(kraken2_trimd_summary),
    path(krona_trimd_file)
    val(extended_qc)
    val(coverage)
    path(nodes_file)
    path(names_file)

    output:
    tuple val(meta), path('*_summary_complete.txt'),            emit: outcome
    path('*_summaryline.tsv'),                                  optional:true, emit: summary_line
    tuple val(meta), path('*.synopsis'),                        optional:true, emit: synopsis
    path("versions.yml"),                                       emit: versions

    script:
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fairy_read_count_outcome_file = fairy_read_count_outcome ? "$fairy_read_count_outcome" : ""
    def raw_qc = raw_qc_file ? "-a $raw_qc_file" : ""
    def fastp_total_qc_pipeline_stats = fastp_total_qc_file ? "-b $fastp_total_qc_file" : ""
    def fastp_total_qc_summaryline = fastp_total_qc_file ? "-t $fastp_total_qc_file" : ""
    def kraken2_trimd_summary_pipeline_stats = kraken2_trimd_summary ? "-f $kraken2_trimd_summary" : ""
    def kraken2_trimd_summary_summaryline = kraken2_trimd_summary ? "-k $kraken2_trimd_summary" : ""
    def kraken2_trimd_report = kraken2_trimd_report_file ? "-e $kraken2_trimd_report_file" : ""
    def krona_trimd = krona_trimd_file ? "-g $krona_trimd_file" : ""
    def extended_qc_arg = extended_qc ? "--extended_qc" : ""
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script_id = params.ica ? "${params.ica_path}/determine_taxID.sh" : "determine_taxID.sh"
    def script_writer = params.ica ? "${params.ica_path}/pipeline_stats_writer.sh" : "pipeline_stats_writer.sh"
    def script_summary = params.ica ? "python ${params.ica_path}/Phoenix_summary_line.py" : "Phoenix_summary_line.py"
    def script_edit = params.ica ? "python ${params.ica_path}/edit_line_summary.py" : "edit_line_summary.py"
    def terra = params.terra ? "-2 terra" : ""
    """
    # set new final script name
    complete_summary="${prefix}_summary_complete.txt"
    
    # handle -entry SCAFFOLDS
    scaffold_entry_file() {
        cat <<<EOT >> \${complete_summary}
        PASSED: Using Scaffold entry no corruption check run on R1.
        PASSED: Using Scaffold entry no corruption check run on R2.
        PASSED: Using Scaffold entry no paired reads to check.
        PASSED: Using Scaffold entry no trimd reads to check.
        FAILED: No scaffolds in ${prefix} after filtering!
        EOT
    }

    # checking that the output contains scaffolds still:
    if grep "Output:                 	0 reads (0.00%) 	0 bases (0.00%)" ${bbmap_log}; then
        #Check if the file exists already (it won't with -entry SCAFFOLDS)
        if [ -f ${fairy_read_count_outcome} ]; then
            # replace end of line with actual error message
            cp ${fairy_read_count_outcome} \${complete_summary}
            sed -i 's/End_of_File/FAILED: No scaffolds in ${prefix} after filtering!/' \${complete_summary}
        else
            scaffold_entry_file
            echo "FAILED: No scaffolds in ${prefix} after filtering!" >> \${complete_summary}
        fi

        # if the sample has no scaffolds left make the summaryline and synopsis file for it. 
        # get taxa ID
        ${script_id} -r $kraken2_trimd_summary -s ${prefix} -d $nodes_file -m $names_file

        # write synopsis file
        ${script_writer} \\ 
            -d ${prefix} \\
            -q ${prefix}.tax \\
            -5 $coverage \\
            $raw_qc \\
            $fastp_total_qc_pipeline_stats \\
            $kraken2_trimd_report \\
            $kraken2_trimd_summary_pipeline_stats \\
            $krona_trimd $terra

        # write summary_line file
        ${script_summary} \\
            -n ${prefix} \\
            -s ${prefix}.synopsis \\
            -x ${prefix}.tax
            -o ${prefix}_summaryline.tsv\\
            $kraken2_trimd_summary_summaryline \\
            $fastp_total_qc_summaryline \\
            $extended_qc_arg

        # change pass to fail and add in error
        ${script_edit} -i ${prefix}_summaryline.tsv

    # if there are scaffolds left after filtering do the following...
    else
        #Check if the file exists already (it won't with -entry SCAFFOLDS)
        if [ -f ${fairy_read_count_outcome} ]; then
            #replace end of line with actual error message
            cp ${fairy_read_count_outcome} \${complete_summary}
            sed -i 's/End_of_File/PASSED: More than 0 scaffolds in ${prefix} after filtering./' \${complete_summary}
        else
            scaffold_entry_file
            echo "PASSED: More than 0 scaffolds in ${prefix} after filtering." >> \${complete_summary}
        fi
    fi

    #gettings script versions
    dettaxid_version=\$(${script_id} -V)
    pipestats_version=\$(${script_writer} -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${dettaxid_version}
        \${pipestats_version}
        Phoenix_summary_line.py: \$(${script_summary} --version )
        edit_line_summary.py: \$(${script_edit} --version )
    END_VERSIONS
    """
}