process GRIPHIN {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(summary_line_files)
    path(original_samplesheet)
    path(db)
    path(outdir)
    val(coverage)

    output:
    path("*_Report.xlsx"),           emit: griphin_report
    path("GRiPHin_samplesheet.csv"), emit: converted_samplesheet
    path("versions.yml"),            emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """
    full_path=\$(readlink -f ${outdir})

    GRiPHin.py -d \$full_path -a $db --output ${outdir} --coverage ${coverage}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}