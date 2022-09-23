def VERSION = '1.0' // Version information not provided by tool on CLI

process GET_MLST_SRST2 {
    tag "${meta.id}"
    label 'process_low'

    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
    //    'quay.io/biocontainers/srst2:0.2.0--py27_2'}"

    //container "quay.io/jvhagey/phoenix:base_v1.0.1"
    container "quay.io/biocontainers/python:2.7--1"

    input:
    tuple val(meta),  path(taxonomy)

    output:
    tuple val(meta), path("*_getMLST_out.txt")                                                , emit: getMLSTs
    tuple val(meta), path("*.fasta")                                           , optional:true, emit: fastas
    tuple val(meta), path("*_profiles.csv")                                    , optional:true, emit: profiles
    tuple val(meta), path("*_pull_dates.txt")                                                 , emit: pull_date
    path "versions.yml"                                                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
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
    python -V
    which python
    echo "here we go"
    convert_MLST_DB_spaces.py --genus "\${genus}" --species "\${species}" > DB_defs.txt
    echo "there we went"

    dbline=\$(tail -n1 DB_defs.txt)
    echo "\$dbline"
    IFS=',' read -r -a db_array <<< "\$dbline"
    echo "\${#db_array[@]}-\${db_array[@]}"
    #num_dbs="\${#db_array[@]}"
    counter=1
    for entry in "\${db_array[@]}"
    do
      echo "Entry#\${counter}-\${entry}|"
      entry_no_spaces="\${entry// /_}"
      if [[ "\${entry}" = "No match found" ]]; then
    		touch "\${entry_no_spaces}.fasta"
    		touch "\${entry_no_spaces}_profiles.csv"
        echo "DB:No match found(\${genus} \${species})       defs:\${entry_no_spaces}_profiles.csv        del:''" > \${entry_no_spaces}_getMLST_out.txt
      else
        getMLST2_phoenix.py --species "\$entry"
        if [[ "\${entry}" = *"baumannii#1" ]]; then
    			sed -i -e 's/Oxf_//g' "\${entry_no_spaces}.fasta"
    			sed -i -e 's/Oxf_//g' "\${entry_no_spaces}_profiles.csv"
    		elif [[ "\${entry}" = *"baumannii#2" ]]; then
    			sed -i -e 's/Pas_//g' "\${entry_no_spaces}.fasta"
    			sed -i -e 's/Pas_//g' "\${entry_no_spaces}_profiles.csv"
        fi
      fi
      echo "\${today}" > "\${entry_no_spaces}_pull_dates.txt"
      counter=\$(( counter + 1))
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        getMLST: $VERSION
    END_VERSIONS
    """
}
