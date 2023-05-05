process GRIPHIN {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(summary_line_files)
    path(original_samplesheet)
    path(db)
    path(outdir)

    output:
    path("*_Report.xlsx"),           emit: griphin_report
    path("GRiPHin_samplesheet.csv"), emit: converted_samplesheet
    path("versions.yml"),            emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """
    #create a samplesheet to be passed to GRiPHin.py
    #while IFS="" read -r line;
    #do
    #    sample_name=\$(echo \$line | cut -d ',' -f 1)
     #   echo 
    #    if [[ "\$sample_name" == "sample" ]]; then
    #        echo "sample,directory" > GRiPHin_samplesheet.csv
    #    else
    #        #get the full path for the samples rather than the working directory
    #        full_path=\$(readlink -f ${outdir}/Phoenix_Output_Report.tsv)
    #        full_dir=\$(echo \$full_path | sed 's/\\/Phoenix_Output_Report.tsv//')
    #        echo \$sample_name,\$full_dir/\$sample_name >> GRiPHin_samplesheet.csv
    #    fi
    #done < ${original_samplesheet}

    full_path=\$(readlink -f ${outdir})

    GRiPHin.py -d \$full_path -a $db --output ${outdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}