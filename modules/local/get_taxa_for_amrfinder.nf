process GET_TAXA_FOR_AMRFINDER {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.0.2'

    input:
    tuple val(meta), path(taxa_file)

    output:
    tuple val(meta), path("*_AMRFinder_Organism.csv"), emit: amrfinder_taxa
    path("versions.yml"),                           emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    get_taxa_for_amrfinder.py -t $taxa_file -o ${prefix}_AMRFinder_Organism.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
