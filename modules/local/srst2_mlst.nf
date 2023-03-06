process SRST2_MLST {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
        'quay.io/biocontainers/srst2:0.2.0--py27_2'}"

    input:
    tuple val(meta), path(fastqs), path(getmlstout), path(alleles), path(profiles), val(status)

    output:
    //tuple val(meta), path("*_mlst_*_results.txt")                , optional:true, emit: mlst_results
    tuple val(meta), path("*_srst2.mlst")                        , optional:true, emit: mlst_results
    tuple val(meta), path("*_srst2_temp.mlst")                   , optional:true, emit: mlst_results_temp
    tuple val(meta), path("*.pileup")                            , optional:true, emit: pileup
    tuple val(meta), path("*.sorted.bam")                        , optional:true, emit: sorted_bam
    tuple val(meta), path("*_srst2_status.txt")                  , optional:true, emit: empty_checker
    path "versions.yml"                                          ,                emit: versions

    when:
    (task.ext.when == null || task.ext.when) //&& "${status[0]}" == "False"

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_s = meta.single_end ? "--input_se ${fastqs}" : "--input_pe ${fastqs[0]} ${fastqs[1]}"
    """
    echo "STATUS-IN: ${status[0]}"
    srst2_ran="False"
    if [[ "${status[0]}" = "False" ]]; then
      scheme_count=1
      for getout in $getmlstout
      do
        no_match="False"
        echo "\${getout}"
        line="\$(tail -n1 \${getout})"
        if [[ "\${line}" = "DB:No match found"* ]] || [[ "\${line}" = "DB:Server down"* ]]; then
          no_match="True"
          mlst_db="No match found"
        else
          # Pulls suggested command info from the getmlst script
          mlst_db=\$(echo "\${line}" | cut -f1 | cut -d':' -f2)
          mlst_defs=\$(echo "\${line}" | cut -f2 | cut -d':' -f2)
          # because this is so messed up and cant pass things through nextflow easily
          #mlst_db="\${mlst_db//.fasta/temp.fasta}"
          #mlst_defs="\${mlst_defs//_profiles.csv/_profiles_temp.csv}"
          mv "\${mlst_db}_temp.fasta" "\${mlst_db}.fasta"
          mv "\${mlst_defs//_profiles.csv/_profiles_temp.csv}" "\${mlst_defs}"
          mlst_delimiter=\$(echo "\${line}" | cut -f3 | cut -d':' -f2 | cut -d"'" -f2)
          #mlst_delimiter=\$(echo "\${line}" | cut -f3 | cut -d':' -f2)

          #mv \${mlst_db}_profiles_csv profiles_csv

          echo "Test: \${mlst_db} \${mlst_defs} \${mlst_delimiter}"

          srst2 ${read_s} \\
              --threads $task.cpus \\
              --output \${scheme_count}_${prefix} \\
              --mlst_db "\${mlst_db}".fasta \\
              --mlst_definitions \${mlst_defs} \\
              --mlst_delimiter \${mlst_delimiter} \\
              $args
        fi
        header="Sample	database	ST	mismatches	uncertainty	depth	maxMAF	locus_1	locus_2	locus_3	locus_4	locus_5	locus_6	locus_7	locus_8	locus_9	locus_10"
        if [[ "\${scheme_count}" -eq 1 ]]; then
          echo "\${header}" > ${prefix}_srst2.mlst
        fi
        header_list=""
        trailer_list=""
        if [[ "\${no_match}" = "True" ]]; then
          tax_with_no_scheme=\$(echo "\${line}" | cut -d'(' -f2 | cut -d')' -f1)
          echo "${prefix}	No match found for \${tax_with_no_scheme}	-	-	-	-	-" >> "${prefix}_srst2.mlst"
        else
          raw_header="\$(head -n1 \${scheme_count}_${prefix}*.txt)"
          raw_trailer="\$(tail -n1 \${scheme_count}_${prefix}*.txt)"
          formatted_trailer="${prefix}	\${mlst_db}"
          IFS=\$'\t' read -r -a trailer_list <<< "\$raw_trailer"
          IFS=\$'\t' read -r -a header_list <<< "\$raw_header"
          header_length="\${#header_list[@]}"
          ST_index=1
          mismatch_index=\$(( header_length - 4 ))
          uncertainty_index=\$(( header_length - 3 ))
          depth_index=\$(( header_length - 2 ))
          maxMAF_index=\$(( header_length - 1))
          genes_start_index=2
          genes_end_index=\$(( header_length - 5 ))
          # Original Index will be something as follows for a 7 gene scheme
          #
          # Sample  ST      adk     atpA    dxr     glyA    recA    sodA    tpi     mismatches      uncertainty     depth   maxMAF
          # 0       1       2       3       4       5       6       7       8       9                10             11      12
          #
          # Expected Index will be as follows for the same isolate to accomomdate splunk ingestion of varying gene counts
          #
          # Sample  database  ST  mismatches  uncertainty depth maxMAF  locus_1 locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 locus_8 locus_9 locus_10
          # 0     1   2   3           4           5     6       7       8       9       10      11      12        13    14      15      16

          echo "\${#header_list[@]} --- \${header_list[@]} --- \${#trailer_list[@]} --- \${trailer_list[@]}"
          formatted_trailer="\${formatted_trailer}	\${trailer_list[\${ST_index}]}"
          formatted_trailer="\${formatted_trailer}	\${trailer_list[\${mismatch_index}]}"
          formatted_trailer="\${formatted_trailer}	\${trailer_list[\${uncertainty_index}]}"
          formatted_trailer="\${formatted_trailer}	\${trailer_list[\${depth_index}]}"
          formatted_trailer="\${formatted_trailer}	\${trailer_list[\${maxMAF_index}]}"

          #for index in {\$genes_start_index..\$genes_end_index}
          for (( index=\${genes_start_index} ; index <= \${genes_end_index} ; index++ ));
          do
            echo "\${index} -- \${header_list[\${index}]} -- \${trailer_list[\${index}]}"
            formatted_trailer="\${formatted_trailer}	\${header_list[\${index}]}(\${trailer_list[\${index}]})"
          done
          echo "\${formatted_trailer}" >> ${prefix}_srst2.mlst
        fi

        #if [[ -f \${scheme_count}_${prefix}*.txt ]]; then
        #  rm \${scheme_count}_${prefix}*.txt
        #fi
        scheme_count=\$(( scheme_count + 1 ))
        #mv profiles_csv \${mlst_db}_profiles_csv
        cp ${prefix}_srst2.mlst ${prefix}_srst2_temp.mlst
        srst2_ran="True"
      done
    else
      #Create a temp file to allow downstream join to complete
      echo "DONT USE" > ${prefix}_srst2_temp.mlst
    fi

    echo "\${srst2_ran}" > "${prefix}_srst2_status.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        srst2: \$(echo \$(srst2 --version 2>&1) | sed 's/srst2 //' )
    END_VERSIONS
    """
}
