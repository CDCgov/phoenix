process KRAKEN2_KRONA {
    tag "$meta.id"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 22)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta), path(kraken_report)
    val(type) //weighted, trimmmed or assembled

    output:
    tuple val(meta), path('*.krona'), emit: krona
    path("versions.yml")            , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/ orginally from https://github.com/jenniferlu717/KrakenTools on 6/15/2022
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.2.0"
    def krakentools_version = "1.2"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}kreport2krona.py \\
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