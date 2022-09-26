process SRST2_MLST {
    tag "${meta.id}"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
        'quay.io/biocontainers/srst2:0.2.0--py27_2'}"

    input:
    tuple val(meta), path(fastqs), path(getmlstout), path(alleles), path(profiles)

    output:
    //tuple val(meta), path("*_mlst_*_results.txt")              , optional:true, emit: mlst_results
    tuple val(meta), path("*_srst2.mlst")                        , optional:true, emit: mlst_results
    tuple val(meta), path("*.pileup")                            , optional:true, emit: pileup
    tuple val(meta), path("*.sorted.bam")                        , optional:true, emit: sorted_bam
    path "versions.yml"                                          ,                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_s = meta.single_end ? "--input_se ${fastqs}" : "--input_pe ${fastqs[0]} ${fastqs[1]}"
    """
    counter=1
    for getout in $getmlstout
    do
      no_match="False"
      echo "\${getout}"
      line="\$(tail -n1 \${getout})"
      if [[ "\${line}" = "DB:No match found"* ]]; then
        no_match="True"
        mlst_db="No match found"
      else
        # Pulls suggested command info from the getmlst script
        mlst_db=\$(echo "\${line}" | cut -f1 | cut -d':' -f2)
        mlst_defs=\$(echo "\${line}" | cut -f2 | cut -d':' -f2)
        mlst_delimiter=\$(echo "\${line}" | cut -f3 | cut -d':' -f2 | cut -d"'" -f2)
        #mlst_delimiter=\$(echo "\${line}" | cut -f3 | cut -d':' -f2)

        #mv \${mlst_db}_profiles_csv profiles_csv

        echo "Test: \${mlst_db} \${mlst_defs} \${mlst_delimiter}"

        srst2 ${read_s} \\
            --threads $task.cpus \\
            --output \${counter}_${prefix} \\
            --mlst_db "\${mlst_db}".fasta \\
            --mlst_definitions \${mlst_defs} \\
            --mlst_delimiter \${mlst_delimiter} \\
            $args
      fi
      header="Sample  database  ST  locus_1 locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 Extra_info(extra_loci,CC,srst2_match_info)"
      if [[ "\${counter}" -eq 1 ]]; then
        if [[ "\${no_match}" = "True" ]];
          tax_with_no_scheme=\$(echo "\${line}" | cut -d'(' -f2 | cut -d')' -f1)
          echo "\${header}" > ${prefix}_srst2.mlst
          echo "${prefix} No match found  - - - - - - - - \${tax_with_no_scheme}" >> "${prefix}_srst2.mlst"
        else
          raw_header="\$(head -n1 1_${prefix}*.txt)"
          raw_trailer="\$(tail -n1 1_${prefix}*.txt)"
          formatted_trailer="${prefix}  \${mlst_db}"
          IFS=' ' read -r -a trailer_list <<< "\$raw_trailer"
          IFS=' ' read -r -a header_list <<< "\$raw_header"
          inner_counter=0
          found_last_locus="False"
          for item in \$trailer_list
          do
            if [[ "\${innner_counter}" -eq 1 ]]; then
              formatted_trailer="\${formatted_trailer}  \${mlst_db}]"
            elif [[ "\${inner_counter}" -ge 2 ]] && [[ "\${inner_counter}" -lt 10 ]]; then
              formatted_trailer="\${formatted_trailer}  \${header_list}[\${inner_counter}](\${trailer_list}[\${inner_counter}])"
            elif [[ "\${inner_counter}" -eq 2 ]] || [[ "\${inner_counter}" -ge 10 ]]; then
              if [[ "\${header_list}[\${inner_counter}]" = "mismatches" ]]; then
                found_last_locus="True"
              fi
              if [[ "\${found_last_locus}" = "False" ]]; then
                formatted_trailer="\${formatted_trailer}  \${header_list}[\${inner_counter}](\${trailer_list}[\${inner_counter}])"
              else
                formatted_trailer="\${formatted_trailer}  \${trailer_list}[\${inner_counter}]"
              fi
            fi
            inner_counter=\$(( inner_counter + 1 ))
          done
          #full_header="database \${header}"
          full_trailer="\${mlst_db} \${formatted_trailer}"
          echo "\${header}" > ${prefix}_srst2.mlst
          echo "\${full_trailer}" >> ${prefix}_srst2.mlst
        fi
      else
        if [[ "\${no_match}" = "True" ]];
          tax_with_no_scheme=\$(echo "\${line}" | cut -d'(' -f2 | cut -d')' -f1)
          echo "${prefix} No match found  - - - - - - - - \${tax_with_no_scheme}" >> "\${counter_${prefix}.txt"
        else
          trailer="\$(tail -n1 \${counter}_${prefix}*.txt)"
          full_trailer="\${mlst_db} \${trailer}"
          echo "\${full_trailer}" >> ${prefix}_srst2.mlst
        fi
      fi
      rm \${counter}_${prefix}*.txt
      counter=\$(( counter + 1 ))
      #mv profiles_csv \${mlst_db}_profiles_csv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        srst2: \$(echo \$(srst2 --version 2>&1) | sed 's/srst2 //' ))
    END_VERSIONS
    """
}
