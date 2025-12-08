process CREATE_NCBI_UPLOAD_SHEET {
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 25)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f7cb3aa4e3324cab43d8635be17da8ae15f62e39d380acda844d1c9deef69c60'

    input:
    path(griphin_samplesheet)
    path(microbe_example)
    path(sra_metadata)
    path(osii_bioprojects)
    path(outdir)
    path(griphin_tsv_report)

    output:
    path("*_Sra_Microbe.1.0.xlsx"),                 optional: true, emit: ncbi_sra
    path("*_BiosampleAttributes_Microbe.1.0.xlsx"), optional: true, emit: ncbi_biosample
    path("versions.yml"),                                           emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}get_ncbi_format_file.py -d ${outdir} --biosample-type microbe -o ./ -s ${sra_metadata} -m ${microbe_example} -b ${osii_bioprojects} -g ${griphin_tsv_report}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        get_ncbi_format_file.py: \$(${ica}get_ncbi_format_file.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}