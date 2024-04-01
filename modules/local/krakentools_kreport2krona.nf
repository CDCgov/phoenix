process KRAKEN2_KRONA {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 22)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(kraken_report)
    val(type) //weighted, trimmmed or assembled

    output:
    tuple val(meta), path('*.krona'), emit: krona
    path("versions.yml")            , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/ orginally from https://github.com/jenniferlu717/KrakenTools on 6/15/2022
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.1.0"
    def krakentools_version = "1.2"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script = params.ica ? "python ${params.ica_path}/kreport2krona.py" : "kreport2krona.py"
    """
    ${script} \\
        --report ${kraken_report} \\
        --output ${prefix}_${type}.krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools: ${krakentools_version}
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}