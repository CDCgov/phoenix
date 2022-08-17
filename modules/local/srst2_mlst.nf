process SRST2_MLST {
    tag "${meta.id}"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/srst2%3A0.2.0--py27_2':
        'quay.io/biocontainers/srst2:0.2.0--py27_2'}"

    input:
    tuple val(meta), path(fastq_s)
    path(taxonomy)

    output:
    tuple val(meta), path("getmlst.out")               , optional:true, emit: gene_results
    tuple val(meta), path("*_mlst_*_results.txt")                , optional:true, emit: mlst_results
    tuple val(meta), path("*.pileup")                            ,                emit: pileup
    tuple val(meta), path("*.sorted.bam")                        ,                emit: sorted_bam
    path "versions.yml"                                          ,                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ""
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_s = meta.single_end ? "--input_se ${fastq_s}" : "--input_pe ${fastq_s[0]} ${fastq_s[1]}"
    """

    species=\$(tail -n2 ${taxonomy} | head -n1 | cut -d\$'\t' -f2)
    genus=\$(tail -n3 ${taxonomy} | head -n1 | cut -d\$'\t' -f2)
    echo "\${genus} ___ \${species}"
    python -V

    db_entry="\${genus} \${species}"

    getMLST2.py --species '\$db_entry' > getmlst.out

    # Pulls suggested command info from the getmlst script
    suggested_command=\$(tail -n2 "getmlst.out" | head -n1)
    mlst_db=\$(echo "\${suggested_command}" | cut -d' ' -f11)
    mlst_defs=\$(echo "\${suggested_command}" | cut -d' ' -f13)
    mlst_delimiter=\$(echo "\${suggested_command}" | cut -d' ' -f15)

    srst2 ${read_s} \\
        --threads $task.cpus \\
        --output ${prefix} \\
        --mlst_db \${mlst_db} \\
        --mlst_definitions \${mlst_defs} \\
        --mlst_delimiter \${mlst_delimiter} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        srst2: \$(echo \$(srst2 --version 2>&1) | sed 's/srst2 //' ))
    END_VERSIONS