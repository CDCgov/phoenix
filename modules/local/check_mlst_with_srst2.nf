def VERSION = '1.1' // Version information not provided by tool on CLI

process CHECK_MLST_WITH_SRST2 {
    tag "$meta.id"
    label 'process_single'
    container "quay.io/jvhagey/phoenix:base_v1.1.0"

    input:
    tuple val(meta), path(mlst_file), path(srst2_file), path(taxonomy_file), val(status)
    path(local_dbases)

    output:
    tuple val(meta), path("*_combined.tsv"), emit: checked_MLSTs
    tuple val(meta), path("*_status.txt"),   emit: status
    path("versions.yml"),                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) {
        terra = ""
    } else if (params.terra==true) {
        terra = "--no-check-certificate"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
    wget $terra --secure-protocol=TLSv1_3 "https://pubmlst.org/data/dbases.xml" -O dbases2.xml

    if [[ ! -s dbases2.xml ]]; then
        mv dbases.xml local_dbases.xml
        mv dbases2.xml dbases.xml
    else
        mv dbases2.xml remote_dbases.xml
    fi

    if [[ "${status[0]}" == "True" ]]; then
        fix_MLST2.py --input $mlst_file --srst2 $srst2_file --taxonomy $taxonomy_file --mlst_database $local_dbases
    elif [[ "${status[0]}" == "False" ]]; then
        fix_MLST2.py --input $mlst_file --taxonomy $taxonomy_file --mlst_database $local_dbases
    else 
        echo "Shouldnt be able to get here, but checking just in case"
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        check_mlst: $VERSION
        python: \$(python --version | sed 's/Python //g')
        pubMLST_db_download_date: \$(date +"%d-%m-%y")
    END_VERSIONS
    """
}
