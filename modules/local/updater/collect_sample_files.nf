process COLLECT_SAMPLE_FILES {
    tag "${meta.id}"
    stageInMode 'copy'
    label 'process_low'
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    tuple val(meta), path(dir)

    output:
    //most things are optional here as checks for specific files needed for an entry point are done in check_directory_samplesheet.py
    tuple val(meta), path("${meta.id}/file_integrity/${meta.id}*_summary.txt"),                    optional: true, emit: fairy_summary
    tuple val(meta), path("${meta.id}/CENTAR/${meta.id}_centar_output.tsv"),                       optional: true, emit: centar_output
    tuple val(meta), path("${meta.id}/ANI/${meta.id}_ShigaPass_summary.csv"),                      optional: true, emit: shigapass_output
    tuple val(meta), path("${meta.id}/fastp_trimd/${meta.id}_1.trim.fastq.gz"),                    optional: true, emit: read1
    tuple val(meta), path("${meta.id}/fastp_trimd/${meta.id}_2.trim.fastq.gz"),                    optional: true, emit: read2
    tuple val(meta), path("${meta.id}/assembly/${meta.id}.filtered.scaffolds.fa.gz"),              optional: true, emit: scaffolds
    tuple val(meta), path("${meta.id}/annotation/${meta.id}.gff"),                                 optional: true, emit: gff
    tuple val(meta), path("${meta.id}/annotation/${meta.id}.faa"),                                 optional: true, emit: faa
    tuple val(meta), path("${meta.id}/gamma_ar/${meta.id}_ResGANNCBI_*_srst2.gamma"),              optional: true, emit: gamma_ar
    tuple val(meta), path("${meta.id}/gamma_hv/${meta.id}_HyperVirulence_*.gamma"),                optional: true, emit: gamma_hv
    tuple val(meta), path("${meta.id}/gamma_pf/${meta.id}_PF-Replicons_*.gamma"),                  optional: true, emit: gamma_pf
    tuple val(meta), path("${meta.id}/AMRFinder/${meta.id}_all_genes{,_*}.tsv"),                   optional: true, emit: amrfinder_report
    tuple val(meta), path("${meta.id}/kraken2_trimd/${meta.id}.kraken2_trimd.top_kraken_hit.txt"), optional: true, emit: kraken_bh
    tuple val(meta), path("${meta.id}/quast/${meta.id}_summary.tsv"),                              optional: true, emit: quast_report
    tuple val(meta), path("${meta.id}/ANI/${meta.id}_REFSEQ_*.fastANI.txt"),                       optional: true, emit: ani_best_hit
    tuple val(meta), path("${meta.id}/ANI/${meta.id}_REFSEQ_*.ani.txt"),                           optional: true, emit: ani
    tuple val(meta), path("${meta.id}/qc_stats/*_trimmed_read_counts.txt"),                        optional: true, emit: trimmed_stats
    tuple val(meta), path("${meta.id}/mlst/*_combined.tsv"),                                       optional: true, emit: combined_mlst
    tuple val(meta), path("${meta.id}/${meta.id}_Assembly_ratio_*.txt"),                           optional: true, emit: assembly_ratio
    tuple val(meta), path("${meta.id}/${meta.id}.synopsis"),                                                       emit: synopsis
    tuple val(meta), path("${meta.id}/${meta.id}.tax"),                                            optional: true, emit: tax
    tuple val(meta), path("${meta.id}/${meta.id}_summaryline.tsv"),                                                emit: summary_line
    tuple val(meta), path("${meta.id}/${meta.id}_updater_log.tsv"),                                optional: true, emit: readme
    //path ("${meta.id}-CSF.csv"),                                                                               emit: collect_file
    path("versions.yml"),                                                                                          emit: versions

    script: 
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    # Nothing happens just getting all of thse files into channels

    #echo "${meta.id},${meta.id}/mlst/${meta.id}_combined.tsv,${meta.id}/file_integrity/${meta.id}_summary.txt,${meta.id}/assembly/${meta.id}.filtered.scaffolds.fa.gz,${meta.id}/${meta.id}.tax" > "${meta.id}-CSF.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}