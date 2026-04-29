process BUSCO {
    tag "$meta.id"
    label 'process_high'
    //6.0.0
    container 'staphb/busco@sha256:2bad2085bebc8e691bf4c51c913902671a56c0e14f2e4750e971e6c3a19ad14a'

    input:
    tuple val(meta), path('tmp_input/*'), path(busco_lineages_path) // path to busco lineages - downloads if not set
    val (lineage)    // Recommended: BUSCO lineages file - downloads if not set
    path(config_file)                      // Optional:    busco configuration file

    output:
    tuple val(meta), path("*-busco.batch_summary.txt")    , emit: batch_summary
    tuple val(meta), path("short_summary.*.txt")          , optional: true, emit: short_summaries_txt
    tuple val(meta), path("short_summary.specific.*.txt") , optional: true, emit: short_summaries_specific_txt
    tuple val(meta), path("short_summary.*.json")         , optional: true, emit: short_summaries_json
    tuple val(meta), path("*-busco")                      , emit: busco_dir
    path "versions.yml"                                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}-${lineage}"
    def busco_lineage = lineage.equals('auto_prok') ? '--auto-lineage-prok' : "--lineage_dataset ${lineage}"
    def busco_lineage_dir = busco_lineages_path ? "--offline --download_path ${busco_lineages_path}" : ''
    def container = task.container.toString() - "quay.io/biocontainers/busco@"
    //set up for terra
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = "PATH=/opt/conda/envs/busco/bin:\$PATH"
        terra_exit = """PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/busco/bin:||')" """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
    #adding python path for running busco on terra
    $terra

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

    busco --cpu $task.cpus --in "\$INPUT_SEQS" --out ${prefix}-busco --force --mode genome  $busco_lineage \\
        $busco_lineage_dir \\

    # clean up
    rm -rf "\$INPUT_SEQS"

    # Move files to avoid staging/publishing issues
    mv ${prefix}-busco/batch_summary.txt ${prefix}-busco.batch_summary.txt
    mv ${prefix}-busco/*/short_summary.*.{json,txt} . || echo "Short summaries were not available: No genes were found."

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$(busco --version 2>&1 | sed 's/^BUSCO //')
        busco_container: ${container}
    END_VERSIONS

    #revert python path back to main envs for running on terra
    $terra_exit
    """
}
