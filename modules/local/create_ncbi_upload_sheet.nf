process CREATE_NCBI_UPLOAD_SHEET {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    path(griphin_samplesheet)
    path(microbe_example)
    path(sra_metadata)
    path(osii_bioprojects)
    path(outdir)

    output:
    path("Sra_Microbe.1.0.xlsx"),                 emit: ncbi_samplesheet
    path("BiosampleAttributes_Microbe.1.0.xlsx"), emit: ncbi_biosample
    path("versions.yml"),                         emit: versions

    script:
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    get_ncbi_format_file.py -d ${outdir} --biosample-type microbe -o ./ -s ${sra_metadata} -m ${microbe_example} -b ${osii_bioprojects}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}