process SRST2_MLST {
    tag "${meta.id}_${getmlst_entry.baseName}"
    label 'process_medium'
    // 0.2.0_patched
    container 'quay.io/jvhagey/srst2@sha256:d4a68baf84c8818b59f334989ccbeea044baf14a79aaf1bd95a1f24f69d0dc5b'

    input:
    tuple val(meta), path(fastqs), path(getmlst_entry), path(alleles), path(profiles), val(status)

    output:
    tuple val(meta), path("*_srst2_line.txt")   , optional:true, emit: mlst_line
    tuple val(meta), path("*.pileup")           , optional:true, emit: pileup
    tuple val(meta), path("*.sorted.bam")       , optional:true, emit: sorted_bam
    path "versions.yml"                                        , emit: versions

    when:
    (task.ext.when == null || task.ext.when)

    script:
    // set up terra variables — unchanged from original
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = """export PYTHONPATH=/opt/conda/envs/srst2/lib/python2.7/site-packages/
        PATH=/opt/conda/envs/srst2/bin:\$PATH
        """
        terra_exit = """export PYTHONPATH=/opt/conda/envs/phoenix/lib/python3.7/site-packages/
        PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/srst2/bin:||')"
        """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    // define variables
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_s = meta.single_end ? "--input_se ${fastqs}" : "--input_pe ${fastqs[0]} ${fastqs[1]}"
    def container = task.container.toString() - "quay.io/jvhagey/srst2@"

    """
    #adding python path for running srst2 on terra
    $terra

    echo "STATUS-IN: ${status[0]}"

    if [[ "${status[0]}" = "False" ]]; then
        no_match="False"
        line="\$(tail -n1 ${getmlst_entry})"

        if [[ "\${line}" = "DB:No match found"* ]] || [[ "\${line}" = "DB:Server down"* ]]; then
            no_match="True"
            mlst_db="No_match_found"
        else
            mlst_db=\$(echo "\${line}" | cut -f1 | cut -d':' -f2)
            mlst_delimiter=\$(echo "\${line}" | cut -f3 | cut -d':' -f2 | cut -d"'" -f2)

            echo "Test: \${mlst_db} \${mlst_db}_profiles.csv \${mlst_delimiter}"

            srst2 ${read_s} \\
                --threads $task.cpus \\
                --output \${mlst_db}_${prefix} \\
                --mlst_db \${mlst_db}_temp.fasta \\
                --mlst_definitions \${mlst_db}_profiles_temp.csv \\
                --mlst_delimiter \${mlst_delimiter} \\
                $args
        fi

        if [[ -f \${mlst_db}_${prefix}__mlst__\${mlst_db}_temp__results.txt ]]; then
            lines_in_result_file=\$(cat \${mlst_db}_${prefix}__mlst__\${mlst_db}_temp__results.txt | wc -l)
        else
            lines_in_result_file=0
            echo "No srst2 result file exists"
        fi

        # --- Produce exactly one output line for THIS scheme (or none, for no_match) ---
        # This mirrors the original loop's 3 cases exactly, minus the header-write concern.
        if [[ "\${no_match}" = "True" ]]; then
            # Original behavior: no_match contributes NO line (matches commented-out original)
            echo "No line written: no_match for this scheme" >&2
        elif [[ "\${lines_in_result_file}" -eq 1 ]]; then
            echo "Not enough was found to even make a guess"
            echo "${prefix}	\${mlst_db}	-	-	-	-	-" > "\${mlst_db}_srst2_line.txt"
        else
            raw_header="\$(head -n1 \${mlst_db}_${prefix}*.txt)"
            to_remove_prefixes=('Pas_' 'Ox_')
            trimmed_header="\${raw_header}"
            for p in \${to_remove_prefixes[@]}; do
                trimmed_header="\${trimmed_header//\${p}/}"
            done
            raw_trailer="\$(tail -n1 \${mlst_db}_${prefix}*.txt)"
            formatted_trailer="${prefix}	\${mlst_db}"
            IFS=\$'\t' read -r -a trailer_list <<< "\$raw_trailer"
            IFS=\$'\t' read -r -a header_list <<< "\$trimmed_header"
            header_length="\${#header_list[@]}"
            ST_index=1
            mismatch_index=\$(( header_length - 4 ))
            uncertainty_index=\$(( header_length - 3 ))
            depth_index=\$(( header_length - 2 ))
            maxMAF_index=\$(( header_length - 1))
            genes_start_index=2
            genes_end_index=\$(( header_length - 5 ))

            formatted_trailer="\${formatted_trailer}	\${trailer_list[\${ST_index}]}"
            formatted_trailer="\${formatted_trailer}	\${trailer_list[\${mismatch_index}]}"
            formatted_trailer="\${formatted_trailer}	\${trailer_list[\${uncertainty_index}]}"
            formatted_trailer="\${formatted_trailer}	\${trailer_list[\${depth_index}]}"
            formatted_trailer="\${formatted_trailer}	\${trailer_list[\${maxMAF_index}]}"

            for (( index=\${genes_start_index} ; index <= \${genes_end_index} ; index++ ));
            do
                formatted_trailer="\${formatted_trailer}	\${header_list[\${index}]}(\${trailer_list[\${index}]})"
            done
            echo "\${formatted_trailer}" > "\${mlst_db}_srst2_line.txt"
        fi
    else
        echo "DONT USE" > "skipped_srst2_line.txt"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        srst2: \$(echo \$(srst2 --version 2>&1) | sed 's/srst2 //' )
        srst2_commit_patched: 73f885f55c748644412ccbaacecf12a771d0cae9
        srst2_container: ${container}
    END_VERSIONS

    #revert python path back to main envs for running on terra
    $terra_exit
    """
}