process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'
    //container 'staphb/ncbi-amrfinderplus:3.10.36'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus%3A3.11.2--h6e70893_0':
        'quay.io/biocontainers/ncbi-amrfinderplus:3.11.2--h6e70893_0' }"

    input:
    tuple val(meta), path(nuc_fasta), val(organism_param), path(pro_fasta), path(gff)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_all_genes.tsv"),                    emit: report
    tuple val(meta), path("${meta.id}_all_mutations.tsv"), optional:true, emit: mutation_report
    path("versions.yml")                                 ,                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ( "${organism_param[0]}" != "No Match Found") {
        organism = "--organism ${organism_param[0]}"
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
        --mutation_all ${prefix}_all_mutations.tsv \\
        $organism \\
        --plus \\
        --database amrfinderdb \\
        --threads $task.cpus > ${prefix}_all_genes.tsv

    sed -i '1s/ /_/g' ${prefix}_all_genes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus_db_version: \$(head amrfinderdb/version.txt)
    END_VERSIONS
    """
}
