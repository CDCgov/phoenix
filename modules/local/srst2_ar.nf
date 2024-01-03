process SRST2_AR {
    tag "${meta.id}"
    label 'process_medium'
    // 0.2.0
    container 'quay.io/jvhagey/srst2@sha256:c53cc041efa44805936233065d3808f08807e8964c3052595a5d93a05784b041'

    input:
    tuple val(meta), path(fastq_s), val(fairy_outcome)
    val(database_type)
    path(db)

    output:
    tuple val(meta), path("*_genes_*_results.txt")                              , emit: gene_results
    tuple val(meta), path("*_fullgenes_*_results.txt")                          , emit: fullgene_results
    tuple val(meta), path("*_mlst_*_results.txt")                , optional:true, emit: mlst_results
    tuple val(meta), path("*.pileup")                            , optional:true, emit: pileup
    tuple val(meta), path("*.sorted.bam")                        , optional:true, emit: sorted_bam
    path("versions.yml")                                         ,                emit: versions

    when:
    //if there are scaffolds left after filtering
    "${fairy_outcome[3]}" == "PASSED: There are reads in ${meta.id} R1/R2 after trimming."

    script:
    // options for running
    if (database_type=="gene") {
        database = "--gene_db ${db}"
    } else if (database_type=="mlst") {
        database = "--mlst_db ${db}"
    } else {
        error "Please set meta.db to either \"gene\" or \"mlst\""
    }
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = """export PYTHONPATH=/opt/conda/envs/srst2/lib/python2.7/site-packages/
        PATH=/opt/conda/envs/srst2/bin:\$PATH
        """
        terra_exit = """export PYTHONPATH=/opt/conda/envs/phoenix/lib/python3.7/site-packages/
        PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/srst2/bin:||')"
        """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    // define variables
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_s = meta.single_end ? "--input_se ${fastq_s}" : "--input_pe ${fastq_s[0]} ${fastq_s[1]}"
    def container = task.container.toString() - "quay.io/jvhagey/srst2@"
    """
    #adding python path for running srst2 on terra
    $terra

    srst2 \\
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
        srst2_commit_with_edits: 73f885f55c748644412ccbaacecf12a771d0cae9
        srst2_container: ${container}
        AMR Combined Database: $db
    END_VERSIONS

    #revert python path back to main envs for running on terra
    $terra_exit
    """
}
