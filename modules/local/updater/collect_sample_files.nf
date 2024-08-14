process COLLECT_SAMPLE_FILES {
    tag "${meta.id}"
    stageInMode 'copy'
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(dir)

    output:
    tuple val(meta), path("${meta.id}/file_integrity/${meta.id}*summary.txt"),                     emit: fairy_summary
    tuple val(meta), path("${meta.id}/fastp_trimd/${meta.id}_1.trim.fastq.gz"),                    emit: read1
    tuple val(meta), path("${meta.id}/fastp_trimd/${meta.id}_2.trim.fastq.gz"),                    emit: read2
    tuple val(meta), path("${meta.id}/assembly/${meta.id}.filtered.scaffolds.fa.gz"),              emit: scaffolds
    tuple val(meta), path("${meta.id}/annotation/${meta.id}.gff"),                                 emit: gff
    tuple val(meta), path("${meta.id}/annotation/${meta.id}.faa"),                                 emit: faa
    tuple val(meta), path("${meta.id}/gamma_hv/${meta.id}_HyperVirulence_*.gamma"),                emit: gamma_hv
    tuple val(meta), path("${meta.id}/gamma_pf/${meta.id}_PF-Replicons_*.gamma"),                  emit: gamma_pf
    tuple val(meta), path("${meta.id}/kraken2_trimd/${meta.id}.kraken2_trimd.top_kraken_hit.txt"), emit: kraken_bh
    tuple val(meta), path("${meta.id}/quast/${meta.id}_summary.tsv"),                              emit: quast_report
    tuple val(meta), path("${meta.id}/ANI/${meta.id}_REFSEQ_*.fastANI.txt"),                       emit: ani_best_hit
    tuple val(meta), path("${meta.id}/qc_stats/*_trimmed_read_counts.txt"),                        emit: trimmed_stats
    tuple val(meta), path("${meta.id}/mlst/*_combined.tsv"),                                       emit: combined_mlst
    tuple val(meta), path("${meta.id}/${meta.id}_Assembly_ratio_*.txt"),                           emit: assembly_ratio
    tuple val(meta), path("${meta.id}/${meta.id}.synopsis"),                                       emit: synopsis
    tuple val(meta), path("${meta.id}/${meta.id}.tax"),                                            emit: tax
    tuple val(meta), path("${meta.id}/${meta.id}_summaryline.tsv"),                                emit: summary_line
    path("versions.yml"),                                                                          emit: versions

    script: 
    // define variables
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    # Nothing happens just getting all of thse files into channels 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}