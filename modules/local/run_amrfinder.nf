process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'
    container 'staphb/ncbi-amrfinderplus:3.10.36'

    /*container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus%3A3.10.23--h17dc2d4_0':
        'quay.io/biocontainers/ncbi-amrfinderplus:3.10.23--h17dc2d4_0' }"*/

    input:
    tuple val(meta), path(fasta), val(organism_param)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_amr_hits.tsv")                     , emit: report
    tuple val(meta), path("${meta.id}_all_mutations.tsv"), optional: true, emit: mutation_report 
    path("versions.yml")                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( "${organism_param[0]}" != "No Match Found") {
        organism = "--organism ${organism_param[0]} --mutation_all ${prefix}_all_mutations.tsv"
    } else {
        organism = ""
    }
    fasta_name = fasta.getName().replace(".gz", "")
    fasta_param = "-n"
    if (meta.containsKey("is_proteins")) {
        if (meta.is_proteins) {
            fasta_param = "-p"
        }
    }
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mkdir amrfinderdb
    tar xzvf $db -C amrfinderdb

    amrfinder \\
        $fasta_param $fasta_name \\
        $organism \\
        $args \\
        --database amrfinderdb \\
        --threads $task.cpus > ${prefix}_amr_hits.tsv

    if [ ! -f ${prefix}_all_mutations.tsv ]; then
        touch ${prefix}_all_mutations.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}
