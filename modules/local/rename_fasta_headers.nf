process RENAME_FASTA_HEADERS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::gamma=2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gamma%3A2.1--hdfd78af_0':
        'quay.io/biocontainers/gamma:2.1--hdfd78af_0' }" // calling this container because it has biopython

    input:
    tuple val(meta), path(assembled_scaffolds)

    output:
    tuple val(meta), path('*.renamed.scaffolds.fa.gz'), emit: renamed_scaffolds
    path "versions.yml"                               , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip $assembled_scaffolds

    rename_fasta_headers.py --input ${prefix}.scaffolds.fa --output ${prefix}.renamed.scaffolds.fa --name ${prefix}

    gzip ${prefix}.renamed.scaffolds.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}