process GRIPHEN {
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(summary_line_files)
    //path(fastp_stats), path(kraken2_trmd_summary), path(kraken2_wtasmbld_summary)
    //path(quast_report), path(ani_file), path(gamma_hv_file), path(gamma_ar_file), path(gamma_pf_file),
    //path(busco_file), path(srst2_ar_file), path(mlst_file), path(assembly_ratio_file)
    path(original_samplesheet)
    path(db)

    output:
    path("OA_Report.xlsx"),            emit: griphen_report
    path("samplesheet_converted.csv"), emit: converted_samplesheet


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

    GRiPhen.py -s samplesheet_converted.csv -a $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}