process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 20)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    path samplesheet

    output:
    path('*.valid.csv'), emit: csv
    path("versions.yml"), emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // define variables
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script = params.ica ? "python ${params.ica_path}/check_samplesheet.py" : "check_samplesheet.py"
    """
    ${script} \\
    $samplesheet \\
    samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        check_samplesheet.py: ${script} --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
