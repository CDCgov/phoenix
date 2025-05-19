process CHECK_SHIGAPASS_TAXA {
    tag "${meta.id}"
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 20)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    tuple val(meta), path(fastani_file), path(ani_file), path(shigapass_file)

    output:
    tuple val(meta), path('edited/*.fastANI.txt'), emit: ani_best_hit
    //tuple val(meta), path("${meta.id}.tax"), emit: tax_file
    path("versions.yml"),                   emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    # when running -entry UPDATE_PHOENIX input will have same name as the output so we will create a directory to store the output
    mkdir -p edited

    #get string to rename file --> Remove "to_check_" from the filename
    new_name=\$(echo "${fastani_file}" | sed 's/to_check_//')

    ${ica}check_taxa.py --format_ani_file ${fastani_file} --shigapass_file ${shigapass_file} --ani_file ${ani_file} --output \${new_name}

    #move output to folder for publishing
    mv \${new_name} edited/\${new_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        check_taxa.py: \$(${ica}check_taxa.py --version)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}