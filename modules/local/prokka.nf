process PROKKA {
    tag "$meta.id"
    label 'process_high'
    container 'staphb/prokka:1.14.5'

    /*container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/prokka:1.14.6--pl526_0' :
        'quay.io/biocontainers/prokka:1.14.6--pl526_0' }"*/

    input:
    tuple val(meta), path(fasta)
    path proteins
    path prodigal_tf

    output:
    tuple val(meta), path("*.gff"), emit: gff //GFF3 format, containing both sequences and annotations
    tuple val(meta), path("*.gbk"), emit: gbk
    tuple val(meta), path("*.fna"), emit: fna //Nucleotide FASTA file of the input contig sequences.
    tuple val(meta), path("*.faa"), emit: faa //Protein FASTA file of the translated CDS sequences.
    tuple val(meta), path("*.ffn"), emit: ffn
    tuple val(meta), path("*.sqn"), emit: sqn
    tuple val(meta), path("*.fsa"), emit: fsa
    tuple val(meta), path("*.tbl"), emit: tbl
    tuple val(meta), path("*.err"), emit: err
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    """
    # Main output unzipped formatted fasta headers lines
    FNAME=\$(basename ${fasta} .gz)
    # Original copy of zipped input fasta
    NFNAME=\$(basename ${fasta} .fa.gz)_original.fa.gz
    # Working copy of unziped input fasta to create and format main prokka input fasta (Cant have cov_x in for downstream)
    NFNAME_U=\$(basename \${NFNAME} .gz)
    mv \${fasta} \${NFNAME}
    gunzip -f \${NFNAME}

    while IFS= read -r line || [ -n "\$line" ]; do
        if [[ "\${line}" = ">"* ]]; then
            no_cov_line=\$(echo "\${line}" | rev | cut -d"_" -f3- | rev)
            echo "\${no_cov_line}" >> "\${FNAME}"
        else
            echo "\${line}" >> "\${FNAME}"
        fi
    done < \${NFNAME_U}

    prokka \\
        $args \\
        --cpus $task.cpus \\
        --prefix $prefix \\
        \$FNAME

    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.gff > ./${prefix}.gff
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.gbk > ./${prefix}.gbk
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.fna > ./${prefix}.fna
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.faa > ./${prefix}.faa
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.ffn > ./${prefix}.ffn
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.sqn > ./${prefix}.sqn
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.fsa > ./${prefix}.fsa
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.tbl > ./${prefix}.tbl
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.err > ./${prefix}.err
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.log > ./${prefix}.log
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.txt > ./${prefix}.txt
    sed 's/NODE/${prefix}/g' ${prefix}/${prefix}.tsv > ./${prefix}.tsv

    #mv ${prefix}/* .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
