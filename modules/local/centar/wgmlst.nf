process WGMLST {
    tag "${meta.id}"
    label 'process_medium'
    container params.wgmlst_container

    input:
    tuple val(meta), path(assembly)
    path(blast_kb)
    val(blastdb)
    val(loci)
    path(qc_kb)

    output:
    tuple val(meta), path("outputs.json"),                                                        emit: outputs
    tuple val(meta), path("stats_calls.json.gz"),                                                 emit: stats
    tuple val(meta), path("allele_calls.bam"),                                                    emit: allele_calls_bam
    tuple val(meta), path("allele_calls.xml.gz"), path("allele_calls.json.gz"),                   emit: allele_calls
    tuple val(meta), path("calls_standard.json.gz"),                                              emit: standard_calls
    tuple val(meta), path("calls_core_standard.csv.gz"), path("calls_core_pcr.csv.gz"),           emit: csv_core
    tuple val(meta), path("calls_accessory_standard.csv.gz"), path("calls_accessory_pcr.csv.gz"), emit: csv_accessory
    path("versions.yml"),                                                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    def container = params.wgmlst_container.toString()
    """
    # Call the real internal scripts to infer the ribotpes
    ngs-run AlleleCalling \
        --sample-id ${meta.id} \
        --publish-dir . \
        --n-threads $task.cpus \
        --assembly ${assembly} \
        --blast-kb.similarity $params.blast_similarity \
        --blast-kb.path ${blast_kb} \
        --blast-kb.db ${blastdb} \
        --blast-kb.loci ${loci} \
        --qc-kb.path ${qc_kb} \
        --organism.genus DEFAULT

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wgMLST_container: ${container}
    END_VERSIONS
    """
}
