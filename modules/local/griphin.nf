process GRIPHIN {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(summary_line_files)
    path(original_samplesheet)
    path(db)

    output:
    path("GRiPhin_Report.xlsx"),       emit: griphin_report
    path("samplesheet_converted.csv"), emit: converted_samplesheet
    path("versions.yml")             , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """
    while IFS="" read -r line;
    do
        sample_name=\$(echo \$line | cut -d ',' -f 1)
        echo 
        if [[ "\$sample_name" == "sample" ]]; then
            echo "sample,directory" > samplesheet_converted.csv
        else
            echo \$sample_name,${params.outdir}/\$sample_name >> samplesheet_converted.csv
        fi
    done < ${original_samplesheet}

    GRiPHin.py -s samplesheet_converted.csv -a $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}