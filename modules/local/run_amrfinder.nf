process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'
    //container 'staphb/ncbi-amrfinderplus:3.10.36'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus%3A3.10.40--h6e70893_1':
        'quay.io/biocontainers/ncbi-amrfinderplus:3.10.40--h6e70893_1' }"

    input:
    tuple val(meta), path(nuc_fasta), val(organism_param), path(pro_fasta), path(gff)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_all_genes.tsv")                       , emit: report
    tuple val(meta), path("${meta.id}_all_mutations.tsv"),   optional: true , emit: mutation_report
    path("versions.yml")                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( "${organism_param[0]}" != "No Match Found") {
        organism = "--organism ${organism_param[0]} --mutation_all ${prefix}_all_mutations.tsv"
    } else {
        organism = ""
    }
    """
    if [[ $nuc_fasta = *.gz ]]; then
        NUC_FNAME=\$(basename ${nuc_fasta} .gz)
        gzip -c -d $nuc_fasta > \$NUC_FNAME
    else
        NUC_FNAME = $nuc_fasta
    fi

    mkdir amrfinderdb
    tar xzvf $db -C amrfinderdb

    amrfinder \\
        --nucleotide \$NUC_FNAME \\
        --protein $pro_fasta \\
        --gff $gff \\
        --annotation_format prokka \\
        $organism \\
        --plus \\
        --database amrfinderdb \\
        --threads $task.cpus > ${prefix}_all_genes.tsv

    sed -i '1s/ /_/g' ${prefix}_all_genes.tsv

    if [ ! -f ${prefix}_all_mutations.tsv ]; then
        touch ${prefix}_all_mutations.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
    END_VERSIONS
    """
}
