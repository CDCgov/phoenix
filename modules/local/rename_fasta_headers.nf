process RENAME_FASTA_HEADERS {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.0.0'

    input:
    tuple val(meta), path(assembled_scaffolds)

    output:
    tuple val(meta), path('*.renamed.scaffolds.fa.gz'), emit: renamed_scaffolds
    path "versions.yml"                               , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    gunzip --force ${assembled_scaffolds}
    unzipped=\$(basename ${assembled_scaffolds} .gz) #adding this in to allow alternative file names with -entry SCAFFOLDS --scaffolds_ext

    rename_fasta_headers.py --input \$unzipped --output ${prefix}.renamed.scaffolds.fa --name ${prefix}

    gzip --force ${prefix}.renamed.scaffolds.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}