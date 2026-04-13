process CREATE_AND_UPDATE_README {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta), path(directory), path(pipeline_info), path(readme), 
      path(old_gamma_ar, stageAs: 'old_gamma_ar/*'), 
      path(new_gamma_ar, stageAs: 'new_gamma_ar/*'), 
      path(old_ncbi_ar,  stageAs: 'old_ncbi_ar/*'), 
      path(new_ncbi_ar,  stageAs: 'new_ncbi_ar/*'), 
      path(old_tax,      stageAs: 'old_tax/*'), 
      path(new_tax,      stageAs: 'new_tax/*'), 
      path(old_pf,       stageAs: 'old_pf/*'), 
      path(new_pf,       stageAs: 'new_pf/*'),
      path(old_software_versions)
    val(current_phx_version)
    path(mlst_db)
    path(ar_db)
    path(pf_db)
    path(amrfinder_db)

    output:
    path('edited/*_updater_log.tsv'), emit: updater_log
    path("versions.yml"),             emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    // allowing for some optional parameters for -entry SCAFFOLDS/CDC_SCAFFOLDS nothing should be passed.
    def old_tax_file = old_tax ? "--old_tax ${old_tax}" : ""
    def new_tax_file = new_tax ? "--new_tax ${new_tax}" : ""
    def old_updater_software_versions = old_software_versions ? "--old_software_version_file ${old_software_versions}" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}Update_Readme.py --mlst_db ${mlst_db} --amrfinder_db ${amrfinder_db} --ar_db ${ar_db} --old_gamma ${old_gamma_ar} --new_gamma ${new_gamma_ar} \\
        --old_ncbi ${old_ncbi_ar} --new_ncbi ${new_ncbi_ar} ${old_tax_file} ${new_tax_file} --old_pf ${old_pf} --new_pf ${new_pf} \\
        -p ${pipeline_info} -d ${directory}/${prefix} -v ${current_phx_version} -o ${prefix}_updater_log.tsv --pf_db ${pf_db} ${old_updater_software_versions}

    #move to output location for process to complete
    mkdir edited/
    mv ${prefix}_updater_log.tsv edited/${prefix}_updater_log.tsv

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
