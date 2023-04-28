def VERSION = '1.2' // Version information not provided by tool on CLI

process KRAKENTOOLS_MAKEKREPORT {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(kraken_output), path(kraken2db_path)

    output:
    tuple val(meta), path('*_wtasmbld.report.txt'), emit: kraken_weighted_report
    path("versions.yml")                          , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/
    // This script has to be run with kraken output that does not use --use-names flag https://github.com/jenniferlu717/KrakenTools/issues/29
    // taxonomy file comes from ./bin/make_ktaxonomy.py --nodes ./taxonomy/nodes.dmp --names ./taxonomy/names.dmp --seqid2taxid seqid2taxid.map --output ktax_map.k2
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    make_kreport.py \\
        --input ${kraken_output} \\
        --output ${prefix}.kraken2_wtasmbld.report.txt \\
        --taxonomy ${kraken2db_path}/ktaxonomy.tsv \\
        --use-read-len
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools_makekreport: $VERSION
    END_VERSIONS
    """
}