process SRATOOLS_PREFETCH {
    tag "${sra_accession[0]}"
    label 'process_single'
    // 3.1.1--h4304569_0
    container "quay.io/biocontainers/sra-tools@sha256:05de2c580cccc4c609ec7c645902563e5d5ffbd366662e1983cb152545ec7bc0"

    input:
    val(sra_accession)

    output:
    path("*_Folder")    , emit: sra_folder
    path('versions.yml'), emit: versions

    script:
    //define variables
    def container = task.container.toString() - "quay.io/biocontainers/sra-tools@"
    """
    # fetch sras
    prefetch --verify yes ${sra_accession[0]} --force all

    #move so we have some common name to collect output, indexing is just to get rid of [] around the SRR number
    mv ${sra_accession[0]} ${sra_accession[0]}_Folder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | sed 's/prefetch : //' | awk 'NF')
        sratools_container: ${container}
    END_VERSIONS
    """
}