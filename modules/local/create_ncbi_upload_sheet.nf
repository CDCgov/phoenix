process CREATE_NCBI_UPLOAD_SHEET {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path(griphin_samplesheet)
    path(microbe_example)
    path(sra_metadata)
    path(osii_bioprojects)
    path(outdir)
    path(griphin_tsv_report)

    output:
    path("Sra_Microbe.1.0.xlsx"),                 emit: ncbi_samplesheet
    path("BiosampleAttributes_Microbe.1.0.xlsx"), emit: ncbi_biosample
    path("versions.yml"),                         emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    ${ica}get_ncbi_format_file.py -d ${outdir} --biosample-type microbe -o ./ -s ${sra_metadata} -m ${microbe_example} -b ${osii_bioprojects} -g ${griphin_tsv_report}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}