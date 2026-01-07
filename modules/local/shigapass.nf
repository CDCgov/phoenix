process SHIGAPASS {
    tag "${meta.id}"
    label 'process_medium'
    container "staphb/shigapass@sha256:83da89164161c54995ec422a55afea39267bc44c194b9a33ccc593ff1d8109e4"

    input:
    tuple val(meta), path(taxa_file), path(scaffolds)
    path(shigapass_database)

    output:
    tuple val(meta), path("*_ShigaPass_Flex_summary.csv"), optional:true, emit: flex_summary
    tuple val(meta), path("*_ShigaPass_summary.csv"),                     emit: summary
    path("versions.yml"),                                                 emit: versions

    script:
    def container = task.container.toString() - "staphb/shigapass@"
    def version = "1.5.0"
    def database_path = shigapass_database ? "${shigapass_database}" : "/ShigaPass-${version}/SCRIPT/ShigaPass_DataBases/"
    """
    #first we create a list of the reads to pass
    gunzip --force ${scaffolds}
    echo ${meta.id}.filtered.scaffolds.fa > ShigaPass_input.txt

    ShigaPass.sh -l ShigaPass_input.txt \\
        -t $task.cpus \\
        -o ShigaPass_Results \\
        -p ${database_path} > Shiga_log.txt

    # catch if the isolate is not Shigella and edit the _ShigaPass_summary.csv so it has the correct header
    if grep -q "Not Shigella/EIEC" Shiga_log.txt; then
        # write header
        echo "Name;rfb;rfb_hits,(%);MLST;fliC;CRISPR;ipaH;Predicted_Serotype;Predicted_FlexSerotype;Comments" > ${meta.id}_ShigaPass_summary.csv 
        echo "${meta.id};;;;;;;;;Not Shigella/EIEC" >> ${meta.id}_ShigaPass_summary.csv
    else
        mv ShigaPass_Results/ShigaPass_summary.csv ./${meta.id}_ShigaPass_summary.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(ShigaPass.sh -v 2>&1 | sed 's/^ShigaPass version //')
        shigapass_container: ${container}
    END_VERSIONS
    """
}