process GRIPHEN {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(fastp_stats), path(kraken2_trmd_summary), path(kraken2_wtasmbld_summary),
    path(quast_report), path(ani_file), path(gamma_hv_file), path(gamma_ar_file), path(gamma_pf_file),
    path(busco_file), path(srst2_ar_file), path(mlst_file), path(assembly_ratio_file)
    path(samplesheet)
    path(db)

    output:
    path("OA_Report.xlsx"),            emit: griphen_report
    path("samplesheet_converted.csv"), emit: converted_samplesheet


    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
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
    done < Caz_Avi_samples.csv

    GRiPhen.py -s samplesheet_converted.csv -a $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}