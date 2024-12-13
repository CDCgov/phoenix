process GRIPHIN {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.0.2'

    input:
    path(summary_line_files)
    path(original_samplesheet)
    path(db)
    path(outdir) // output directory used as prefix for the summary file
    val(coverage)
    val(entry)
    val(scaffolds_entry)
    val(update)
    val(shigapass_detected)

    output:
    path("*_GRiPHin*.xlsx"),                 emit: griphin_report
    path("*_GRiPHin*.tsv"),                  emit: griphin_tsv_report
    path("Directory_samplesheet.csv"),       emit: converted_samplesheet
    path("versions.yml"),                    emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def phoenix = entry ? "--phoenix" : ""
    def scaffolds = scaffolds_entry ? "--scaffolds" : ""
    def shigapass = shigapass_detected ? "--shigapass" : ""
    def output_prefix = update ? "${outdir}_GRiPHin" : "${outdir}_GRiPHin_Summary"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    full_path=\$(readlink -f ${outdir})

    ${ica}GRiPHin.py -d \$full_path -a $db --output ${output_prefix} --coverage ${coverage} ${phoenix} ${shigapass} ${scaffolds} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
       griphin.py: \$(${ica}GRiPHin.py --version)
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}