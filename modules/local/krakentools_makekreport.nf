process KRAKENTOOLS_MAKEKREPORT {
    tag "$meta.id"
    label 'process_single'
    // base_v2.3.0 - MUST manually change below (line 22)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(kraken_output), path(kraken2db_path)

    output:
    tuple val(meta), path('*_wtasmbld.summary.txt'), emit: kraken_weighted_report
    path("versions.yml")                          , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/
    // This script has to be run with kraken output that does not use --use-names flag https://github.com/jenniferlu717/KrakenTools/issues/29
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.3.0"
    def krakentools_version = "1.2"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}make_kreport.py \\
        --input ${kraken_output} \\
        --output ${prefix}.kraken2_wtasmbld.summary.txt \\
        --taxonomy ${kraken2db_path}/ktaxonomy.tsv \\
        --use-read-len
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools: ${krakentools_version}
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}