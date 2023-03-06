def VERSION = '1.1' // Version information not provided by tool on CLI

process CHECK_MLST {
    tag "$meta.id"
    label 'process_low'
    container "quay.io/jvhagey/phoenix:base_v1.1.0"

    input:
    tuple val(meta), path(mlst_file), path(taxonomy_file)

    output:
    tuple val(meta), path("*_combined.tsv")                                                   , emit: checked_MLSTs
    tuple val(meta), path("*_status.txt")                                                     , emit: status
    path "versions.yml"                                                                       , emit: versions

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
    wget $terra --secure-protocol=TLSv1_3 "https://pubmlst.org/data/dbases.xml"

    check_and_fix_MLST2_new2.py --input $mlst_file --taxonomy $taxonomy_file --docfile dbases.xml
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        check_mlst: $VERSION
        pubMLST_db_download_date:
    END_VERSIONS
    """
}
