process KRAKENTOOLS_KREPORT2MPA {
    tag "$meta.id"
    label 'process_single'
    // base_v2.3.0 - MUST manually change below (line 21)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(kraken_report)

    output:
    tuple val(meta), path('*.mpa'), emit: mpa
    path("versions.yml")          , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.3.0"
    def krakentools_version = "1.2"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}kreport2mpa.py \\
        --report-file ${kraken_report} \\
        --output ${prefix}.mpa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools: ${krakentools_version}
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}