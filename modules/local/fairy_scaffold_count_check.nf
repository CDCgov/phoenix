process SCAFFOLD_COUNT_CHECK {
    tag "${meta.id}"
    label 'process_medium'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(bbmap_log), path(fairy_read_count_outcome),
    path(raw_qc),
    path(fastp_total_qc),
    path(kraken2_trimd_report),
    path(kraken2_trimd_summary),
    path(krona_trimd)
    val(extended_qc)
    val(coverage)
    path(nodes_file)
    path(names_file)

    output:
    tuple val(meta), path('*_summary.txt'),                emit: outcome
    path('*_summaryline.tsv'),      optional:true, emit: summary_line
    tuple val(meta), path('*.synopsis'),    optional:true, emit: synopsis
    path("versions.yml"),                                  emit: versions

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
    def extended_qc_arg = extended_qc ? "--extended_qc" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    #checking that the output contains scaffolds still:
    if grep "Output:                 	0 reads (0.00%) 	0 bases (0.00%)" ${bbmap_log}; then
        #echo "FAILED: No scaffolds in ${prefix} after filtering!" >> ${fairy_read_count_outcome}
        #replace end of line with actual error message
        sed -i 's/End_of_File/FAILED: No scaffolds in ${prefix} after filtering!/' ${fairy_read_count_outcome}

        # get taxa ID
        ${ica_bash}determine_taxID.sh -r $kraken2_trimd_summary -s ${prefix} -d $nodes_file -m $names_file

        #write synopsis file
        ${ica_bash}pipeline_stats_writer.sh \\
        -a $raw_qc \\
        -b $fastp_total_qc \\
        -d ${prefix} \\
        -e $kraken2_trimd_report \\
        -f $kraken2_trimd_summary \\
        -g $krona_trimd \\
        -q ${prefix}.tax \\
        -5 $coverage \\
        $terra

        # write summary_line file
        ${ica_python}Phoenix_summary_line.py \\
        -n ${prefix} \\
        -k $kraken2_trimd_summary \\
        -t $fastp_total_qc \\
        -s ${prefix}.synopsis \\
        -x ${prefix}.tax \\
        -o ${prefix}_summaryline.tsv \\
        $extended_qc_arg

        #change file name. 
        cp ${fairy_read_count_outcome} ${prefix}_summary.txt
    else
        sed -i 's/End_of_File/PASSED: More than 0 scaffolds in ${prefix} after filtering./' ${fairy_read_count_outcome}
        #echo "PASSED: More than 0 scaffolds in ${prefix} after filtering." >> ${fairy_read_count_outcome}
        cp ${fairy_read_count_outcome} ${prefix}_summary.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
