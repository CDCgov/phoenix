process GET_TAXA_FOR_AMRFINDER {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/python:3.8.3' :
      'quay.io/biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(taxa_file)

    output:
    tuple val(meta), path("*_AMRFinder_Organism.csv"), emit: amrfinder_taxa

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    get_taxa_for_amrfinder.py -t $taxa_file -o ${prefix}_AMRFinder_Organism.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
