process BUSCO {
    tag "$meta.id"
    label 'process_high'
    // 5.4.7--pyhdfd78af_0
    container 'quay.io/biocontainers/busco@sha256:f5ef1f64076deb66ed015d0b6692619610eb1de22f0a9741bbf0ea8434d06404'

    input:
    tuple val(meta), path('tmp_input/*'), path(busco_lineages_path), // path to busco lineages - downloads if not set
    val(fairy_outcome)
    each(lineage)                          // Required:    lineage to check against, "auto" enables --auto-lineage instead
    path(config_file)                      // Optional:    busco configuration file

    output:
    tuple val(meta), path("*-busco.batch_summary.txt")    , emit: batch_summary
    tuple val(meta), path("short_summary.*.txt")          , optional: true, emit: short_summaries_txt
    tuple val(meta), path("short_summary.specific.*.txt") , optional: true, emit: short_summaries_specific_txt
    tuple val(meta), path("short_summary.*.json")         , optional: true, emit: short_summaries_json
    tuple val(meta), path("*-busco")                      , emit: busco_dir
    path "versions.yml"                                   , emit: versions

    when:
    //if there are scaffolds left after filtering
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def busco_config = config_file ? "--config $config_file" : ''
    def busco_lineage = lineage.equals('auto') ? '--auto-lineage' : "--lineage_dataset ${lineage}"
    def busco_lineage_dir = busco_lineages_path ? "--offline --download_path ${busco_lineages_path}" : ''
    def container = task.container.toString() - "quay.io/biocontainers/busco@"
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = "export PYTHONPATH=/opt/conda/envs/busco/lib/python3.7/site-packages/"
        terra_exit = "export PYTHONPATH=/opt/conda/envs/phoenix/lib/python3.7/site-packages/"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
    #adding python path for running busco on terra
    $terra

    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi

    # If the augustus config directory is not writable, then copy to writeable area
    if [ ! -w "\${AUGUSTUS_CONFIG_PATH}" ]; then
        # Create writable tmp directory for augustus
        AUG_CONF_DIR=\$( mktemp -d -p \$PWD )
        cp -r \$AUGUSTUS_CONFIG_PATH/* \$AUG_CONF_DIR
        export AUGUSTUS_CONFIG_PATH=\$AUG_CONF_DIR
        echo "New AUGUSTUS_CONFIG_PATH=\${AUGUSTUS_CONFIG_PATH}"
    fi

    # Ensure the input is uncompressed
    INPUT_SEQS=input_seqs
    mkdir "\$INPUT_SEQS"
    cd "\$INPUT_SEQS"
    for FASTA in ../tmp_input/*; do
        if [ "\${FASTA##*.}" == 'gz' ]; then
            gzip -cdf "\$FASTA" > \$( basename "\$FASTA" .gz )
        else
            ln -s "\$FASTA" .
        fi
    done
    cd ..

    busco \\
        --cpu $task.cpus \\
        --in "\$INPUT_SEQS" \\
        --out ${prefix}-busco \\
        $busco_lineage \\
        $busco_lineage_dir \\
        $busco_config \\
        $args \\

    # clean up
    rm -rf "\$INPUT_SEQS"

    # Move files to avoid staging/publishing issues
    mv ${prefix}-busco/batch_summary.txt ${prefix}-busco.batch_summary.txt
    mv ${prefix}-busco/*/short_summary.*.{json,txt} . || echo "Short summaries were not available: No genes were found."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
        busco_container: ${container}
    END_VERSIONS

    #revert python path back to main envs for running on terra
    $terra_exit
    """
}
