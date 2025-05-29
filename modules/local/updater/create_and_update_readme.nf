process CREATE_AND_UPDATE_README {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    tuple val(meta), path(directory), path(pipeline_info)
    val(current_phx_version)
    path(mlst_db)
    path(ar_db)
    path(amrfinder_db)

    output:
    path('*_updater_log.tsv'), emit: updater_log
    path("versions.yml"),      emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    // allowing for some optional parameters for -entry SCAFFOLDS/CDC_SCAFFOLDS nothing should be passed.
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Update_Readme.py \\
        --mlst_db ${mlst_db} \\
        --amrfinder_db ${amrfinder_db} \\
        --ar_db ${ar_db} \\
        -p ${pipeline_info} \\
        -d ${directory}/${prefix} \\
        -v ${current_phx_version} \\
        -o ${prefix}_updater_log.tsv

    # since Readme is just being appended we need to create a file to signal we are done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        Update_Readme.py: \$(${ica}Update_Readme.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
