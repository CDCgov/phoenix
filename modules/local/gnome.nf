process GNOME {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(combined_raw_stats)
    tuple val(meta), path(reads)

    output:
    tuple val(meta), val('*.fastq.gz'), emit: bbduk_reads

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [[gnome.py ${combined_raw_stats} == PASS]]:
    then
        ${reads}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}