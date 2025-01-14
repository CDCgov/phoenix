process CHECK_SHIGAPASS_TAXA {
    tag "${meta.id}"
    label 'process_low'
    // base_v2.1.0 - MUST manually change below (line 27)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(fastani_file), path(ani_file), path(shigapass_file)

    output:
    tuple val(meta), path('*.fastANI.txt'), emit: ani_best_hit
    //tuple val(meta), path("${meta.id}.tax"), emit: tax_file
    path("versions.yml"),                   emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    #get string to rename file --> Remove "to_check_" from the filename
    new_name=\$(echo "${fastani_file}" | sed 's/to_check_//')

    ${ica}check_taxa.py --format_ani_file ${fastani_file} --shigapass_file ${shigapass_file} --ani_file ${ani_file} --output \${new_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        check_taxa.py: \$(${ica}check_taxa.py --version)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}