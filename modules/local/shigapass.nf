process SHIGAPASS {
    tag "${meta.id}"
    label 'process_medium'
    container "quay.io/jvhagey/newtype@sha256:a002c7e7aa137b31fbd244831d0a3c583f3c2d55b831ef805242bbce44ae97b7"

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
    def container = task.container.toString() - "jvhagey/newtype:"
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

    # catch if the isolate is not Shigella
    if grep -q "Not Shigella/EIEC" Shiga_log.txt; then
        # write header
        echo "Name;rfb;rfb_hits,(%);MLST;fliC;CRISPR;ipaH;Predicted_Serotype;Predicted_FlexSerotype;Comments" > ${meta.id}_ShigaPass_summary.csv 
        echo "${meta.id};;;;;;;;;Not Shigella/EIEC" >> ${meta.id}_ShigaPass_summary.csv
    else
        mv ShigaPass_Results/ShigaPass_summary.csv ./${meta.id}_ShigaPass_summary.csv
        #get taxa
        if grep -iq "sonnei" filename.txt; then
            species="s:624\tsonnei"
        elif grep -iq "flexneri" filename.txt; then
            species="s:623\tflexneri"
        fi
        printf "ShigaPass\t \t${meta.id}_ShigaPass_summary.csv\n" > ${meta.id}.tax
        printf "K:2\tBacteria\nP:1224\tPseudomonadota\nC:1236\tGammaproteobacteria\nO:91347\tEnterobacterales\nF:543\tEnterobacteriaceae\nG:620\tShigella\n\${species}" >> ${meta.id}.tax
        if [ -s ShigaPass_Results/ShigaPass_Flex_summary.csv ]; then
            mv ShigaPass_Results/ShigaPass_Flex_summary.csv ./${meta.id}_ShigaPass_Flex_summary.csv #only makes for S. flexneri genomes
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^ShigaPass version //' ))
        shigapass_container: ${container}
    END_VERSIONS
    """
}