def VERSION = 'https://github.com/jenniferlu717/KrakenTools/commit/fb1fcc4b747470f33b798f6a602d736b9f78c881' // Version information not provided by tool on CLI

process KRAKENTOOLS_MAKEKREPORT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(kraken_output)
    path(ktax_map)

    output:
    tuple val(meta), path('*_wtasmbld.report.txt'), emit: kraken_weighted_report
    path "versions.yml"                           , emit: versions

    script: // This script is bundled with the pipeline, in phoenix/bin/
    // This script has to be run with kraken output that does not use --use-names flag https://github.com/jenniferlu717/KrakenTools/issues/29
    // see  folder /scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/minikraken2DB/
    // taxonomy file comes from ./QuAISAR_Nextflow/bin/make_ktaxonomy.py --nodes ./taxonomy/nodes.dmp --names ./taxonomy/names.dmp --seqid2taxid seqid2taxid.map --output ktax_map.k2
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    make_kreport.py \\
        --input ${kraken_output} \\
        --output ${prefix}.kraken2_wtasmbld.report.txt \\
        --taxonomy ${ktax_map} \\
        --use-read-len
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        krakentools_makekreport: $VERSION
    END_VERSIONS
    """
}