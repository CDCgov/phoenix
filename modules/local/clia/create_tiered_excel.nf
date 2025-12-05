process CREATE_TIERED_EXCEL {
    tag "$meta.id"
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(amrfinder_files), path(abritamr_files)

    output:
    tuple val(meta), path("*_AR_details.csv"), emit: tiered_report
    path("versions.yml"),                      emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    #make csv file
    ${ica}process_amr_data.py

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        process_amr_data.py: \$(${ica}process_amr_data.py --version)
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}