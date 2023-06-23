process SRST2_AR {
    tag "${meta.id}"
    label 'process_medium'
    //container 'staphb/srst2:0.2.0'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
        'quay.io/biocontainers/srst2:0.2.0--py27_2'}"

    input:
    tuple val(meta), path(fastq_s), path(db)

    output:
    tuple val(meta), path("*_genes_*_results.txt")                              , emit: gene_results
    tuple val(meta), path("*_fullgenes_*_results.txt")                          , emit: fullgene_results
    tuple val(meta), path("*_mlst_*_results.txt")                , optional:true, emit: mlst_results
    tuple val(meta), path("*.pileup")                            , optional:true, emit: pileup
    tuple val(meta), path("*.sorted.bam")                        , optional:true, emit: sorted_bam
    path "versions.yml"                                          ,                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_s = meta.single_end ? "--input_se ${fastq_s}" : "--input_pe ${fastq_s[0]} ${fastq_s[1]}"
    if (meta.db=="gene") {
        database = "--gene_db ${db}"
    } else if (meta.db=="mlst") {
        database = "--mlst_db ${db}"
    } else {
        error "Please set meta.db to either \"gene\" or \"mlst\""
    }
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = "export PYTHONPATH=/opt/conda/envs/srst2/lib/python2.7/site-packages/"
        terra_exit = "export PYTHONPATH=/opt/conda/envs/phoenix/lib/python3.7/site-packages/"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
    #adding python path for running srst2 on terra
    $terra

    srst2.py \\
        ${read_s} \\
        --threads $task.cpus \\
        --output ${prefix} \\
        ${database} \\
        $args

    # create an empty fullgenes file if nothing was found, otherwise the pre-summary join fails silently
    short_DB=\$(basename ${db} .fasta)
    if [[ ! -f ${prefix}__fullgenes__\${short_DB}__results.txt ]]; then
        touch ${prefix}__fullgenes__\${short_DB}__results.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        srst2: \$(echo \$(srst2 --version 2>&1) | sed 's/srst2 //' )
        AMR Combined Database: $db
    END_VERSIONS

    #revert python path back to main envs for running on terra
    $terra_exit
    """
}
