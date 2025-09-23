process CLIA_GRIPHIN {
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 28)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    path(db)
    path(original_samplesheet)
    val(metas)
    path(griphin_files) // path to the GRiPHin files, which includes the staged directory
    path(outdir) // output directory used as prefix for the summary file
    val(coverage)

    output:
    path("*_GRiPHin_Summary.xlsx"),    emit: griphin_report
    path("*_GRiPHin_Summary.tsv"),     emit: griphin_tsv_report
    path("Phoenix_Summary.tsv"),       emit: phoenix_tsv_report
    path("Directory_samplesheet.csv"), emit: converted_samplesheet
    path("versions.yml"),              emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def shigapass = shigapass_detected ? "--shigapass" : ""
    def prefix = task.ext.prefix ?: "GRiPHin"
    def stage_files = [
        metas.collect { "mkdir -p ${prefix}/${it.id}" },
        metas.collect { "mv ${it.filenames.join(' ')} ${prefix}/${it.id}" }
    ].flatten().join(" && ")
    """
    full_path=\$(readlink -f ${outdir})

    ${ica}CLIA_GRiPHin.py -d \$full_path -a $db --output ${outdir} --coverage ${coverage} 
            --phx_version ${phx_version} ${shigapass}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        griphin.py: \$(${ica}GRiPHin.py --version)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}