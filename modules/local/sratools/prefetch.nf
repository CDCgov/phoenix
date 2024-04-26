process SRATOOLS_PREFETCH {
    tag "${sra_accession[0]}"
    label 'process_single'
    // 3.0.3--h87f3376_0 "quay.io/biocontainers/sra-tools@sha256:c9f92683e10091c3ef93066b1fcbdeeba89af49242ab778a9c8cc006f6be82a3"
    // 3.0.9--h9f5acd7_0
    // 3.1.0--h4304569_1
    container "quay.io/biocontainers/sra-tools@sha256:8c0e26d2ee39a93860168327113674127297ce0cd9bf33ebd489a4fdef08f112"

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
    prefetch --verify yes ${sra_accession[0]}

    #move so we have some common name to collect output, indexing is just to get rid of [] around the SRR number
    mv ${sra_accession[0]} ${sra_accession[0]}_Folder

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(prefetch --version 2>&1 | sed 's/prefetch : //' | awk 'NF')
        sratools_container: ${container}
    END_VERSIONS
    """
}