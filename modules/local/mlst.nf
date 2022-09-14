process MLST {
    tag "$meta.id"
    label 'process_low'
    container 'staphb/mlst:2.22.1' // version 2.22.1 produces add 

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) {
        terra = ""
    } else if (params.terra==true) {
        terra = "PATH=/opt/conda/envs/phoenix/bin/wget:$PATH"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
    $terra
    mlst \\
        --threads $task.cpus \\
        $fasta \\
        > ${prefix}.tsv

    scheme=\$(cut -d \$'\t' -f2 ${prefix}.tsv)
    if [ \$scheme == "abaumannii_2" ]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        sed -i 's/abaumannii_2/abaumannii(Pasteur)/' ${prefix}_1.tsv
        mlst --scheme abaumannii --threads $task.cpus $fasta > ${prefix}_2.tsv
        sed -i 's/abaumannii/abaumannii(Oxford)/' ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [ \$scheme == "abaumannii" ]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        sed -i 's/abaumannii/abaumannii(Oxford)/' ${prefix}_1.tsv
        mlst --scheme abaumannii_2 --threads $task.cpus $fasta > ${prefix}_2.tsv
        sed -i 's/abaumannii/abaumannii(Pasteur)/' ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [ \$scheme == "ecoli_achtman_4" ]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        sed -i 's/ecoli_achtman_4/ecoli(Achtman)/' ${prefix}_1.tsv
        mlst --scheme ecoli --threads $task.cpus $fasta > ${prefix}_2.tsv
        sed -i 's/ecoli/ecoli(Pasteur)/' ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [ \$scheme == "ecoli" ]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        sed -i 's/ecoli/ecoli(Pasteur)/' ${prefix}_1.tsv
        mlst --scheme ecoli_achtman_4 --threads $task.cpus $fasta > ${prefix}_2.tsv
        sed -i 's/ecoli_achtman_4/ecoli(Achtman)/' ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    else
        :
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
    """

}