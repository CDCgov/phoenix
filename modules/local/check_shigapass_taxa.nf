process CHECK_SHIGAPASS_TAXA {
    tag "${meta.id}"
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 20)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(fastani_file), path(ani_file), path(shigapass_file), path(tax_file)

    output:
    tuple val(meta), path('edited/*.fastANI.txt'),              emit: ani_best_hit
    tuple val(meta), path('edited/*.ani.txt'),                  emit: ani_raw
    tuple val(meta), path("edited/${meta.id}.tax"),             emit: tax_file
    tuple val(meta), path("${meta.id}_updater_log.tax"),        emit: edited_tax_file
    path("versions.yml"),                                       emit: versions

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    # when running --mode UPDATE_PHOENIX input will have same name as the output so we will create a directory to store the output
    mkdir -p edited

    # In regular Phoenix mode, fastani_file is FORMAT_ANI's "to_check_" copy, and stripping
    # that prefix yields the name of a genuine, separate sibling file that already exists
    # on disk -- check_taxa.py writes corrections there, leaving fastani_file untouched.
    #
    # In UPDATE_PHOENIX mode, FORMAT_ANI never runs, so there is no "to_check_" prefix to
    # strip and no sibling file waiting to receive corrections -- new_name would otherwise
    # resolve to the SAME path as fastani_file, causing check_taxa.py to overwrite the only
    # existing copy in place. To keep this module correct in both modes without needing to
    # detect which one is active, corrections are always written to an explicitly distinct
    # temp filename instead of relying on new_name being different from fastani_file.
    new_name=\$(echo "${fastani_file}" | sed 's/to_check_//')
    python_output="corrected_\${new_name}"

    ${ica}check_taxa.py --format_ani_file ${fastani_file} --shigapass_file ${shigapass_file} --ani_file ${ani_file} --format_ani_output \${python_output} --tax_file ${tax_file}

    # check_taxa.py still computes and writes corrected values to \${python_output} (used to
    # derive the percent-identity provenance recorded in the .tax file), but we intentionally
    # publish the ORIGINAL fastani_file and the ORIGINAL, untouched ani_file here -- the .tax
    # file is the authoritative final call and records where it came from; the ANI-derived
    # files are left as a historical record of what the raw comparison actually showed.
    cp ${fastani_file} edited/\${new_name}
    cp ${ani_file} edited/
    cp ${meta.id}.tax edited/${meta.id}.tax
    cp edited/${meta.id}.tax ${meta.id}_updater_log.tax # renaming so there isn't a file name conflict when we create the updater log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        \$(${ica}check_taxa.py --version)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}