process CDIFF_RIBOTYPER {
    tag "$meta.id"
    label 'process_single'
    // phx_ml_v1.0.0 - MUST manually change below (line 27)!!!
    container 'quay.io/jvhagey/phoenix@sha256:827455416027d49d1a5f37bdea1acd2fc5c4ab537fa59a97b7ea31df334e48cc'

    input:
    tuple val(meta), path(csv_core)
    tuple val(meta), path(csv_accessory)
    path(newtype_bin_dir)

    output:
    tuple val(meta), path("*_ribotype.tsv"), emit: ribotype_file
    path("versions.yml")                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    def container_version = "phx_ml_v1.0.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    # Call the real internal scripts to infer the ribotpes
    ${ica}newtype.py -i ./ -s PN2.0 -o ${prefix}_ribotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        newtype.py: \$(${ica}${ica}newtype.py --version)
        phoenix_ml_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
