process SHIGAPASS {
    tag "${meta.id}"
    label 'process_medium'
    container "staphb/shigapass@sha256:83da89164161c54995ec422a55afea39267bc44c194b9a33ccc593ff1d8109e4"

    input:
    tuple val(meta), path(scaffolds), path(taxa_file), val(fairy_outcome)
    path(shigapass_database)

    output:
    tuple val(meta), path("*_ShigaPass_Flex_summary.csv"), optional:true, emit: flex_summary
    tuple val(meta), path("*_ShigaPass_summary.csv"),                     emit: summary
    tuple val(meta), path("${meta.id}.tax"),                              emit: tax_file
    path("versions.yml"),                                                 emit: versions

    when:
    //if there are scaffolds left after filtering
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

    script:
    def container = task.container.toString() - "staphb/shigapass@"
    def version = "1.5.0"
    def database_path = shigapass_database ? "-p ${shigapass_database}" : "/ShigaPass-${version}/SCRIPT/ShigaPass_DataBases/"
    """
    #first we create a list of the reads to pass
    gunzip --force ${scaffolds}
    echo ${meta.id}.filtered.scaffolds.fa > ShigaPass_input.txt

    ShigaPass.sh -l ShigaPass_input.txt \\
        -t $task.cpus \\
        -o ShigaPass_Results \\
        -p ${database_path} > Shiga_log.txt

    ## Set a default value for species to avoid unbound variable errors 
    species="s:unknown\tunknown"

    # catch if the isolate is not Shigella
    if grep -q "Not Shigella/EIEC" Shiga_log.txt; then
        # write header
        echo "Name;rfb;rfb_hits,(%);MLST;fliC;CRISPR;ipaH;Predicted_Serotype;Predicted_FlexSerotype;Comments" > ${meta.id}_ShigaPass_summary.csv 
        echo "${meta.id};;;;;;;;;Not Shigella/EIEC" >> ${meta.id}_ShigaPass_summary.csv
    else
        mv ShigaPass_Results/ShigaPass_summary.csv ./${meta.id}_ShigaPass_summary.csv
        #get taxa
        if grep -q "SS" ${meta.id}_ShigaPass_summary.csv; then
            species="s:624\tsonnei"
        elif grep -q "SF1-5" ${meta.id}_ShigaPass_summary.csv; then
            species="s:623\tflexneri"
        elif grep -q "SB" ${meta.id}_ShigaPass_summary.csv; then
            species="s:621\tboydii"
        elif grep -q "SD" ${meta.id}_ShigaPass_summary.csv; then
            species="s:622\tdysenteriae"
        fi
        echo "ShigaPass\t\t${meta.id}_ShigaPass_summary.csv" > ${meta.id}.tax
        echo "K:2\tBacteria\\nP:1224\tPseudomonadota\\nC:1236\tGammaproteobacteria\\nO:91347\tEnterobacterales\\nF:543\tEnterobacteriaceae\\nG:620\tShigella\\n\${species}" >> ${meta.id}.tax
        if [ -s ShigaPass_Results/ShigaPass_Flex_summary.csv ]; then
            mv ShigaPass_Results/ShigaPass_Flex_summary.csv ./${meta.id}_ShigaPass_Flex_summary.csv #only makes for S. flexneri genomes
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(ShigaPass.sh -v 2>&1 | sed 's/^ShigaPass version //')
        shigapass_container: ${container}
    END_VERSIONS
    """
}