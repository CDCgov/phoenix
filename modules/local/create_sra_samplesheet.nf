process CREATE_SRA_SAMPLESHEET {
    label 'process_single'
    // base_v2.3.0 - MUST manually change below (line 23)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    path(renamed_reads)
    path(metadata_csvs)
    path(directory)
    val(srr_param)

    output:
    path('sra_samplesheet.csv'), emit: csv
    path("versions.yml"),        emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def use_srr = srr_param ? "--use_srr" : ""
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    full_path=\$(readlink -f ${directory})

    ${ica}sra_samplesheet.py -d \$full_path $use_srr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sra_samplesheet.py: \$(${ica}sra_samplesheet.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}