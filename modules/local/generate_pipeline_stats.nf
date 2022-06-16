process GENERATE_PIPELINE_STATS {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(fastp_raw_qc)
    tuple val(meta), path(fastp_total_qc)
    tuple val(meta), path(trimmed_reads)
    tuple val(meta), path(kraken2_trimd_report)
    tuple val(meta), path(kraken2_trimd_summary)
    tuple val(meta), path(krona_trimd)
    tuple val(meta), path(assembly_scaffolds)
    tuple val(meta), path(filtered_assembly)
    tuple val(meta), path(kraken2_asmbld_report)
    tuple val(meta), path(kraken2_asmbled_summary)
    tuple val(meta), path(krona_asmbld)
    tuple val(meta), path(kraken2_weighted_report)
    tuple val(meta), path(kraken2_weighted_summary)
    tuple val(meta), path(krona_weighted)
    tuple val(meta), path(quast_report)
    tuple val(meta), path(taxID)
    tuple val(meta), path(assembly_ratio_file)
    tuple val(meta), path(busco_specific_short_summary)
    tuple val(meta), path(fastANI_formatted_file)
    tuple val(meta), path(gamma_AR)
    tuple val(meta), path(gamma_replicon)
    tuple val(meta), path(gamma_HV)
    tuple val(meta), path(srst_fullgenes_file)
    tuple val(meta), path(mlst_file)

    output:
    tuple val(meta), path('*.synopsis'), emit: pipeline_stats

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    pipeline_stats_writer.sh \\
        -a $fastp_raw_qc \\
        -b $fastp_total_qc \\
        -c ${trimmed_reads[0]} \\
        -d ${trimmed_reads[1]} \\
        -e $kraken2_trimd_report \\
        -f $kraken2_trimd_summary \\
        -g $krona_trimd\\
        -h $assembly_scaffolds \\
        -i $filtered_assembly \\
        -j $kraken2_asmbld_report \\
        -k $kraken2_asmbled_summary \\
        -l $krona_asmbld \\
        -m $kraken2_weighted_report \\
        -n $kraken2_weighted_summary \\
        -o $krona_weighted \\
        -p $quast_report \\
        -q $taxID \\
        -r $assembly_ratio_file \\
        -s ${busco_specific_short_summary[1]} \\
        -t $fastANI_formatted_file \\
        -u $gamma_AR \\
        -v $gamma_replicon \\
        -w $gamma_HV \\
        -x $srst_fullgenes_file \\
        -y $mlst_file
    """
}