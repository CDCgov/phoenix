process CREATE_TIERED_EXCEL {
    tag "$meta.id"
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

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