process CREATE_SRA_SAMPLESHEET {
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    path(renamed_reads)
    path(metadata_csvs)
    path(directory)
    val(srr_param)

    output:
    path('sra_samplesheet.csv'), emit: csv
    path("versions.yml"),        emit: versions

    script:
    def use_srr = srr_param ? "--use_srr" : ""
    """
    full_path=\$(readlink -f ${directory})

    sra_samplesheet.py -d \$full_path $use_srr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}