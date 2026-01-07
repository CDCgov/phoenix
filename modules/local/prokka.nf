process PROKKA {
    tag "$meta.id"
    label 'process_high'
    // 1.14.5
    container 'staphb/prokka@sha256:4ef8e13b87f6ba1bc79f599970ec25c60a80913ab0bc15c90171c9743e86994f'

    input:
    tuple val(meta), path(fasta)
    path(proteins)
    path(prodigal_tf)

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

    script:
    //set up for terra
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = "PATH=/opt/conda/envs/prokka/bin:\$PATH"
        terra_exit = """PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/prokka/bin:||')" """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    //define variables
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def proteins_opt = proteins ? "--proteins ${proteins[0]}" : ""
    def prodigal_opt = prodigal_tf ? "--prodigaltf ${prodigal_tf[0]}" : ""
    def container = task.container.toString() - "staphb/prokka@"
    """
    #adding prokka path for running prokka on terra
    $terra

    # Main output unzipped formatted fasta headers lines
    FNAME=\$(basename ${fasta} .gz)
    # Original copy of zipped input fasta
    NFNAME=\$(basename ${fasta} .fa.gz)_original.fa.gz
    # Working copy of unzipped input fasta to create and format main prokka input fasta (Cant have cov_x in for downstream)
    NFNAME_U=\$(basename \${NFNAME} .gz)
    mv ${fasta} \${NFNAME}
    gunzip -f \${NFNAME}
    
    sed 's/${prefix}/NODE/g' \${NFNAME_U} > \${FNAME}

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
        prokka_container: ${container}
    END_VERSIONS

    #revert path back to main envs for running on terra
    $terra_exit
    """
}