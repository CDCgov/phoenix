process CREATE_CLIA_PDF {
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 23)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    path(directory)
    path(phoenix_summary)
    path(griphin_tsv_summary)
    val(start_time)
    path(ar_database)
    val(coverage)

    output:
    path("WGS_Run_Summary_report*.pdf"),  emit: pdf_summary
    path("WGS_Run_Summary_report*.html"), emit: html_summary
    path("versions.yml"),                 emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/griphin/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}report_gen.py -p ${directory} --phx_version v2.2.0 -t ${start_time} --ar_database ${ar_database} --amrfinder_version 3.12.8 --coverage ${coverage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        report_gen.py: \$(${ica}report_gen.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}