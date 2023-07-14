process SRATOOLS_PREFETCH {
    tag "${sra_accession[0]}"
    label 'process_single'
    container "https://depot.galaxyproject.org/singularity/sra-tools%3A3.0.3--h87f3376_0"

    input:
    val(sra_accession)

    output:
    path("*_Folder")    , emit: sra_folder
    path('versions.yml'), emit: versions

    script:
    """
    # fetch sras
    prefetch --verify yes ${sra_accession[0]}

    #move so we have some common name to collect output, indexing is just to get rid of [] around the SRR number
    mv ${sra_accession[0]} ${sra_accession[0]}_Folder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | sed 's/prefetch : //' | awk 'NF')
    END_VERSIONS
    """
}