process SCAFFOLD_COUNT_CHECK {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

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
    tuple val(meta), path('*_summary.txt'),             emit: outcome
    path('*_summaryline.tsv'),           optional:true, emit: summary_line
    tuple val(meta), path('*.synopsis'), optional:true, emit: synopsis
    path("versions.yml"),                               emit: versions

    script:
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) { terra = ""} 
    else if (params.terra==true) { terra = "-2 terra" }
    else { error "Please set params.terra to either \"true\" or \"false\"" }
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { 
        ica_python = "" 
        ica_bash = ""
    } else if (params.ica==true) { 
        ica_python = "python ${workflow.launchDir}/bin/" 
        ica_bash = "bash ${workflow.launchDir}/bin/" 
    }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
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
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    #checking that the output contains scaffolds still:
    if grep "Output:                 	0 reads (0.00%) 	0 bases (0.00%)" ${bbmap_log}; then
        #Check if the file exists already (it won't with -entry SCAFFOLDS)
        if [ -f ${prefix}_summary_old_3.txt ]; then
            #replace end of line with actual error message
            sed -i 's/End_of_File/FAILED: No scaffolds in ${prefix} after filtering!/' ${fairy_read_count_outcome_file}
        else
            echo "PASSED: Using Scaffold entry no corruption check run on R1." > ${prefix}_summary_old_3.txt
            echo "PASSED: Using Scaffold entry no corruption check run on R2." >> ${prefix}_summary_old_3.txt
            echo "PASSED: Using Scaffold entry no paired reads to check." >> ${prefix}_summary_old_3.txt
            echo "PASSED: Using Scaffold entry no trimd reads to check." >> ${prefix}_summary_old_3.txt
            echo "FAILED: No scaffolds in ${prefix} after filtering!" >> ${prefix}_summary_old_3.txt
        fi

        # if the sample has no scaffolds left make the summaryline and synopsis file for it. 
        # get taxa ID
        ${ica_bash}determine_taxID.sh -r $kraken2_trimd_summary -s ${prefix} -d $nodes_file -m $names_file

        #write synopsis file
        ${ica_bash}pipeline_stats_writer.sh -d ${prefix} -q ${prefix}.tax -5 $coverage $raw_qc $fastp_total_qc_pipeline_stats \\
        $kraken2_trimd_report $kraken2_trimd_summary_pipeline_stats $krona_trimd $terra

        # write summary_line file
        ${ica_python}Phoenix_summary_line.py -n ${prefix} -s ${prefix}.synopsis -x ${prefix}.tax -o ${prefix}_summaryline.tsv\\
        $kraken2_trimd_summary_summaryline $fastp_total_qc_summaryline $extended_qc_arg

        # change pass to fail and add in error
        ${ica_python}edit_line_summary.py -i ${prefix}_summaryline.tsv

        #change file name.
        cp ${prefix}_summary_old_3.txt ${prefix}_summary.txt

    # if there are scaffolds left after filtering do the following...
    else
        #Check if the file exists already (it won't with -entry SCAFFOLDS)
        if [ -f ${prefix}_summary_old_3.txt ]; then
            #replace end of line with actual error message
            sed -i 's/End_of_File/PASSED: More than 0 scaffolds in ${prefix} after filtering./' ${fairy_read_count_outcome_file}
        else
            echo "PASSED: Using Scaffold entry no corruption check run on R1." > ${prefix}_summary_old_3.txt
            echo "PASSED: Using Scaffold entry no corruption check run on R2." >> ${prefix}_summary_old_3.txt
            echo "PASSED: Using Scaffold entry no paired reads to check." >> ${prefix}_summary_old_3.txt
            echo "PASSED: Using Scaffold entry no trimd reads to check." >> ${prefix}_summary_old_3.txt
            echo "PASSED: More than 0 scaffolds in ${prefix} after filtering." >> ${prefix}_summary_old_3.txt
        fi
        cp ${prefix}_summary_old_3.txt ${prefix}_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
