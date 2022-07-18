process MLST {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::mlst=2.19.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mlst:2.19.0--hdfd78af_1' :
        'quay.io/biocontainers/mlst:2.19.0--hdfd78af_1' }"

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
    """
    mlst \\
        --threads $task.cpus \\
        $fasta \\
        > ${prefix}.tsv

    scheme=\$(cut -d \$'\t' -f2 ${prefix}.tsv)
    if [ \$scheme == "abaumannii_2" ]; then
        sed -i 's/abaumannii_2/abaumannii_2(Pasteur)/' ${prefix}.tsv > ${prefix}.tsv
        mlst --scheme abaumannii --threads $task.cpus $fasta >> ${prefix}.tsv
        sed -i 's/abaumannii/abaumannii(Oxford)/' ${prefix}.tsv >> ${prefix}.tsv
    elif [ \$scheme == "abaumannii" ]; then
        sed -i 's/abaumannii/abaumannii(Oxford)/' ${prefix}.tsv > ${prefix}.tsv
        mlst --scheme abaumannii_2 --threads $task.cpus $fasta >> ${prefix}.tsv
        sed -i 's/abaumannii_2/abaumannii_2(Pasteur)/' ${prefix}.tsv >> ${prefix}.tsv
    elif [ \$scheme == "ecoli" ]; then
        sed -i 's/ecoli/ecoli(Achtman)/' ${prefix}.tsv > ${prefix}.tsv
        mlst --scheme ecoli_2 --threads $task.cpus $fasta >> ${prefix}.tsv
        sed -i 's/ecoli_2/ecoli_2(Pasteur)/' ${prefix}.tsv >> ${prefix}.tsv
    elif [ \$scheme == "ecoli_2" ]; then
        sed -i 's/ecoli_2/ecoli_2(Pasteur)/' ${prefix}.tsv > ${prefix}.tsv
        mlst --scheme ecoli --threads $task.cpus $fasta >> ${prefix}.tsv
        sed -i 's/ecoli/ecoli(Achtman)/' ${prefix}.tsv >> ${prefix}.tsv
    else
        :
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
    END_VERSIONS
    """

}