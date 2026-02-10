process GET_MLST_SRST2 {
    tag "${meta.id}"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta),  path(taxonomy), val(status), path(local_mlst_db)

    output:
    //tuple val(meta), path("*_getMLST_out.txt")                                 , optional:true, emit: getMLSTs
    //tuple val(meta), path("*.fasta")                                           , optional:true, emit: fastas
    //tuple val(meta), path("*_profiles.csv")                                    , optional:true, emit: profiles
    tuple val(meta), path("*_getMLST_out_temp.txt")                            , optional:true, emit: getMLSTs_checker
    tuple val(meta), path("*_temp.fasta")                                      , optional:true, emit: fastas_checker
    tuple val(meta), path("*_profiles_temp.csv")                               , optional:true, emit: profiles_checker
    tuple val(meta), path("*_pull_dates.txt")                                  , optional:true, emit: pull_date
    path "versions.yml"                                                                       , emit: versions

    when:
    (task.ext.when == null || task.ext.when) //& "${status[0]}" == "False"

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = task.container.toString() - "quay.io/biocontainers/python@"
    """
    if [[ "${status[0]}" == "False" ]]; then
        genus="empty"
        species="empty"
        today=\$(date '+%Y-%m-%d')
        test_title=\$(tail -n2 ${taxonomy} | head -n1)
        echo "-\${test_title}-"
        if [[ \$test_title = "G:"* ]]; then
            species=\$(tail -n1 ${taxonomy} | cut -f2)
            genus=\$(tail -n2 ${taxonomy} | head -n1 | cut -f2)
        elif [[ \$test_title = "s:"* ]]; then
            species=\$(tail -n2 ${taxonomy} | head -n1 | cut -f2)
            genus=\$(tail -n3 ${taxonomy} | head -n1 | cut -f2)
        else
            echo "-\$test_title-"
        fi
        echo "\${genus}___\${species}"
        # Old way, now use provided DB with different name format
        # convert_taxonomy_with_complexes_to_pubMLST.py --genus "\${genus}" --species "\${species}" > DB_defs.txt
        ${ica}local_MLST_converter.py --genus "\${genus}" --species "\${species}" > DB_defs.txt

        dbline=\$(tail -n1 DB_defs.txt)
        echo "\$dbline"
        IFS=',' read -r -a db_array <<< "\$dbline"
        echo "\${#db_array[@]}-\${db_array[@]}"
        #num_dbs="\${#db_array[@]}"
        counter=1
        for DBID in "\${db_array[@]}"
        do
            echo "Entry#\${counter}-\${DBID}|"
            if [[ "\${DBID}" = "No match found" ]]; then
                touch "No_match_found.fasta"
                touch "No_match_found_profiles.csv"
                touch "No_match_found_pull_dates.txt"
                touch "No_match_found_temp.fasta"
                touch "No_match_found_profiles_temp.csv"
                echo "DB:No match found(\${genus} \${species})       defs:No_match_found_profiles.csv        del:''" > No_match_found_getMLST_out.txt
                cp "No_match_found_getMLST_out.txt" "No_match_found_getMLST_out_temp.txt"
                new_pull_date=\$(head "${local_mlst_db}/db_version")
                echo "\${new_pull_date}" >> "No_match_found_pull_date.txt"
            else
                # Have not found any other delimiters other than underscore in all DB folders
                echo -e "DB:\${DBID}\tdefs:\${DBID}_profile.csv\tdel:'_'" > "\${DBID}_getMLST_out.txt"
                cp "\${DBID}_getMLST_out.txt" "\${DBID}_getMLST_out_temp.txt"
                
                for allele_file in "${local_mlst_db}/pubmlst/\${DBID}/"*".tfa"; do
                    echo "\${allele_file}"
                    cat "\${allele_file}" >> "\${DBID}_temp.fasta"
                done

                cp "${local_mlst_db}/pubmlst/\${DBID}/\${DBID}.txt" "\${DBID}_profiles_temp.csv"
                new_pull_date=\$(head "${local_mlst_db}/db_version")
                echo "\${new_pull_date}" >> "\${DBID}_pull_date.txt"

                if [[ "\${DBID}" = "abaumannii" ]]; then
                    sed -i -e 's/Oxf_//g' "\${DBID}_temp.fasta"
                    sed -i -e 's/Oxf_//g' "\${DBID}_profiles_temp.csv"
                elif [[ "\$DBID}" = "abaumannii_2" ]]; then
                    sed -i -e 's/Pas_//g' "\${DBID}.fasta"
                    sed -i -e 's/Pas_//g' "\${DBID}_profiles.csv"
                fi
            fi
            counter=\$(( counter + 1))
        done
    else
        echo "Did not run" > "${prefix}_getMLST_out_temp.txt"
        echo "Did not run" > "${prefix}_temp.fasta"
        echo "Did not run" > "${prefix}_profiles_temp.csv"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        local_MLST_converter.py: \$(${ica}local_MLST_converter.py --version)
        python: \$(python --version | sed 's/Python //g')
        python_container: ${container_version}
    END_VERSIONS
    """
}
