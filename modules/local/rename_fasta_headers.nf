process RENAME_FASTA_HEADERS {
    tag "$meta.id"
    label 'process_low'
    container 'staphb/gamma:2.1'

    input:
    tuple val(meta), path(assembled_scaffolds)

    output:
    tuple val(meta), path('*.renamed.scaffolds.fa.gz'), emit: renamed_scaffolds
    path "versions.yml"                               , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip --force $assembled_scaffolds

    rename_fasta_headers.py --input ${prefix}.scaffolds.fa --output ${prefix}.renamed.scaffolds.fa --name ${prefix}

    gzip --force ${prefix}.renamed.scaffolds.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}