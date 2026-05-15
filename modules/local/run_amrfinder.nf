process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'
    if (params.mode_upper == "CLIA") { // in case they diverge in the future. 
        //4.2.7-2026-03-24.1
        container 'staphb/ncbi-amrfinderplus@sha256:90e57a65bde22d3270ca13d43e59470aa2aa544b0b377fad7d8b2d1032e9741f'
    } else {
        //4.2.7-2026-03-24.1
        container 'staphb/ncbi-amrfinderplus@sha256:90e57a65bde22d3270ca13d43e59470aa2aa544b0b377fad7d8b2d1032e9741f'
    }

    input:
    tuple val(meta), path(nuc_fasta), val(organism_param), path(pro_fasta), path(gff)
    path(db)

    output:
    tuple val(meta), path("${meta.id}_all_genes_*.tsv"),                    emit: report
    tuple val(meta), path("${meta.id}_all_mutations_*.tsv"), optional:true, emit: mutation_report
    path("versions.yml"),                                                   emit: versions

    script:
    // use --organism
    if ( "${organism_param[0]}" != "No Match Found") {
        organism = "--organism ${organism_param[0]}"
    } else { organism = "" }
    //set up for terra
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = "PATH=/opt/conda/envs/amrfinderplus/bin:\$PATH"
        terra_exit = """PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/amrfinderplus/bin:||')" """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    // define variables
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "staphb/ncbi-amrfinderplus@"
    //get name of amrfinder database file
    def db_name = db.toString() - '.tar.gz'
    def db_date = db_name.find(/\d{8}/)
    """
    #adding python path for running srst2 on terra
    $terra

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
        --mutation_all ${prefix}_all_mutations_${db_date}.tsv \\
        $organism \\
        --plus \\
        --database $db_name \\
        --threads $task.cpus > ${prefix}_all_genes_${db_date}.tsv

    sed -i '1s/ /_/g' ${prefix}_all_genes_${db_date}.tsv

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
