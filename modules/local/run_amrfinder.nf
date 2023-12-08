process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'
    // 3.11.11-2023-04-17.1@sha256:194eec0c758f92c3c8a8884b9f1ddbfb7626977459e8938a6ece98aceb8e3bbd
    // 3.11.26-2023-11-15.1 - new
    container 'staphb/ncbi-amrfinderplus@sha256:2311400ff8576f04e77a250e0665daa05f11b67f6a5901c8eac267378dccc932'

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
    // use --organism
    if ( "${organism_param[0]}" != "No Match Found") {
        organism = "--organism ${organism_param[0]}"
    } else { organism = "" }
    if (params.terra==false) {
        terra_activate = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra_activate = "micromamba activate amrfinderplus"
        terra_exit = "micromamba deactivate"
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    // define variables
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "staphb/ncbi-amrfinderplus@"
    //get name of amrfinder database file
    db_name = db.toString() - '.tar.gz'
    """
    #adding python path for running srst2 on terra
    $terra_activate

    if [[ $nuc_fasta = *.gz ]]; then
        NUC_FNAME=\$(basename ${nuc_fasta} .gz)
        gzip -c -d $nuc_fasta > \$NUC_FNAME
    else
        NUC_FNAME = $nuc_fasta
    fi

    # decompress the amrfinder database
    tar xzvf $db

    amrfinder \\
        --nucleotide \$NUC_FNAME \\
        --protein $pro_fasta \\
        --gff $gff \\
        --annotation_format prokka \\
        --mutation_all ${prefix}_all_mutations.tsv \\
        $organism \\
        --plus \\
        --database $db_name \\
        --threads $task.cpus > ${prefix}_all_genes.tsv

    sed -i '1s/ /_/g' ${prefix}_all_genes.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus_db_version: \$(head $db_name/version.txt)
        amrfinderplus_container: ${container} 
    END_VERSIONS

    #revert python path back to main envs for running on terra
    $terra_exit
    """
}
