process CHECK_SHIGAPASS_TAXA {
    tag "${meta.id}"
    label 'process_low'
    // base_v2.3.0 - MUST manually change below (line 20)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(fastani_file), path(ani_file), path(shigapass_file), path(tax_file)

    output:
    tuple val(meta), path('edited/*.fastANI.txt'),              emit: ani_best_hit
    tuple val(meta), path("edited/${meta.id}.tax"),             emit: tax_file
    tuple val(meta), path("${meta.id}_updater_log.tax"),        emit: edited_tax_file
    path("versions.yml"),                                       emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    # when running --mode UPDATE_PHOENIX input will have same name as the output so we will create a directory to store the output
    mkdir -p edited

    #get string to rename file --> Remove "to_check_" from the filename
    new_name=\$(echo "${fastani_file}" | sed 's/to_check_//')

    # Check if the shigella species in the shigapass file matches the species in the fastani file
    if grep -q "Shigella" "${tax_file}"; then
        # Extract species from s: line
        fastani_species=\$(grep "^s:" "${tax_file}" | cut -f2)
        # Get second line from summary (split by semicolon)
        shigapass_species=\$(echo "${shigapass_file}" | sed -n '2p' | cut -d';' -f10)
        # this should catch cases of 
        if [[ "\$shigapass_species" != *"\$fastani_species"* ]]; then
            echo "Shigapass species: \$shigapass_species and FastANI species: \$fastani_species do NOT match. Updating taxa file."
            ${ica}check_taxa.py --format_ani_file ${fastani_file} --shigapass_file ${shigapass_file} --ani_file ${ani_file} --format_ani_output \${new_name} --tax_file ${tax_file}

            #After updating files move output to folder for publishing
            mv \${new_name} edited/\${new_name}
            mv ${meta.id}.tax edited/${meta.id}.tax
            cp edited/${meta.id}.tax ${meta.id}_updater_log.tax # renaming so there isn't a file name conflict when we create the updater log
        else 
            echo "Shigella found, and Shigapass species: \$shigapass_species and FastANI species: \$fastani_species match."
            #No changes to taxa files needed just move to output to folder for publishing
            mv ${fastani_file} edited/\${new_name}
            mv ${meta.id}.tax edited/${meta.id}.tax 
            cp edited/${meta.id}.tax ${meta.id}_updater_log.tax # renaming so there isn't a file name conflict when we create the updater log
        fi
    # If fastani said Escherichia AND shigapass did not say "Not Shigella/EIEC" or EIEC we need to update the taxa file --> not sure this would ever happen
    elif grep -q "Escherichia" "${tax_file}" && ! grep -q "EIEC" "${shigapass_file}"; then
        # Extract genera from s: line
        fastani_genera=\$(grep "^G:" "${tax_file}" | cut -f2)
        # Get second line from summary (split by semicolon)
        shigapass_org=\$(echo "${shigapass_file}" | sed -n '2p' | cut -d';' -f10)
        echo "Escherichia found, and Shigapass taxa: \$shigapass_org and FastANI genera: \$fastani_genera do not match. Updating taxa file."
        ${ica}check_taxa.py --format_ani_file ${fastani_file} --shigapass_file ${shigapass_file} --ani_file ${ani_file} --format_ani_output \${new_name} --tax_file ${tax_file}

        # Move output to folder for publishing
        mv \${new_name} edited/\${new_name}
        mv ${meta.id}.tax edited/${meta.id}.tax
        cp edited/${meta.id}.tax ${meta.id}_updater_log.tax # renaming so there isn't a file name conflict when we create the updater log
    else
        echo "Escherichia or Shigella were not found, PHoeNIx filters are broken please open a github ticket https://github.com/CDCgov/phoenix/issues."
        # Should actually not ever get here since we filter so only Escherichia and Shigella enter SHIGAPASS module, but just in case
        mv ${fastani_file} edited/\${new_name}
        mv ${meta.id}.tax edited/${meta.id}.tax
        cp edited/${meta.id}.tax ${meta.id}_updater_log.tax # renaming so there isn't a file name conflict when we create the updater log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        check_taxa.py: \$(${ica}check_taxa.py --version)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}