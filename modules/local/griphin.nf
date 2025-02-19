process GRIPHIN {
    label 'process_low'
    //container 'quay.io/jvhagey/phoenix:base_v2.0.2'
    // base_v2.1.0 - MUST manually change below (line 22)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    path(summary_line_files)
    path(original_samplesheet)
    path(db)
    path(outdir) // output directory used as prefix for the summary file
    val(coverage)
    val(entry)
    val(scaffolds_entry)
    val(updater_entry)
    val(run_centar)
    val(long_read_entry)

    output:
    tuple path("full_path_file.txt"), path("*_GRiPHin*.xlsx"), emit: griphin_report
    path("*_GRiPHin*.tsv"),                                    emit: griphin_tsv_report
    path("Directory_samplesheet.csv"),                         emit: converted_samplesheet
    path("versions.yml"),                                      emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def phoenix = entry ? "--phoenix" : ""
    def scaffolds = scaffolds_entry ? "--scaffolds" : ""
    def centar = run_centar ? "--centar" : ""
    def long_read = long_read_entry ? "--long_read" : ""
    def output_prefix = updater_entry ? "${outdir}_GRiPHin" : "${outdir}_GRiPHin_Summary"
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    full_path=\$(readlink -f ${outdir})

    # Save the value of full_path to a file (this file will be captured in the output block)
    echo \$full_path > full_path_file.txt

    ${ica}GRiPHin.py -d \$full_path -a $db --output ${output_prefix} --coverage ${coverage} ${phoenix} ${centar} ${scaffolds} ${long_read}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
       griphin.py: \$(${ica}GRiPHin.py --version)
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}