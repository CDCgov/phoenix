process SRATOOLS_PREFETCH {
    tag "${sra_accession[0]}"
    label 'process_single'
    // 3.2.0--h4304569_0
    container "quay.io/biocontainers/sra-tools@sha256:db636fa5785c482fe69d836b2c4e24c9a912b9557ed069cad3835d4234b9354e"

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