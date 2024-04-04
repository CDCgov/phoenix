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

    output:
    path("*_Summary.xlsx"),            emit: griphin_report
    path("*_Summary.tsv"),             emit: griphin_tsv_report
    path("Directory_samplesheet.csv"), emit: converted_samplesheet
    path("versions.yml"),              emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // define variables
    def phoenix = entry ? "--phoenix" : ""
    def scaffolds = scaffolds_entry ? "--scaffolds" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    def script = params.ica ? "python ${params.ica_path}/GRiPHin.py" : "GRiPHin.py"
    """
    full_path=\$(readlink -f ${outdir})

    ${script} -d \$full_path -a $db --output ${outdir} --coverage ${coverage} ${phoenix} ${scaffolds}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
       griphin.py: \$(${script} --version)
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}