process SHIGAPASS {
    tag "${meta.id}"
    label 'process_medium'
    container "quay.io/jvhagey/newtype:1.5.0"

    input:
    tuple val(meta), path(scaffolds), val(fairy_outcome)

    output:
    tuple val(meta), path("*_ShigaPass_Flex_summary.csv"), emit: flex_summary
    tuple val(meta), path("*_ShigaPass_summary.csv"),      emit: summary
    path "versions.yml"                             ,      emit: versions

    when:
    //if there are scaffolds left after filtering
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

    script:
    def container = task.container.toString() - "jvhagey/newtype:"
    """
    #first we create a list of the reads to pass
    gunzip ${scaffolds}
    echo ${meta.id}_filtered.scaffolds.fa > ShigaPass.txt

    ShigaPass.sh -l ShigaPass.txt \\
     -t $task.cpus \\
     -o ShigaPass_Results -p ShigaPass_DataBases -u

    mv ShigaPass_Results/ShigaPass_summary.csv ./${meta.id}_ShigaPass_summary.csv
    mv ShigaPass_Results/ShigaPass_Flex_summary.csv ./${meta.id}_ShigaPass_Flex_summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shigapass: \$(echo \$(ShigaPass.sh -v 2>&1) | sed 's/^ShigaPass version //' ))
        shigapass_container: ${container}
    END_VERSIONS
    """
}