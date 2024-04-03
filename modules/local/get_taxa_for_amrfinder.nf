process GET_TAXA_FOR_AMRFINDER {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 21)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(taxa_file)

    output:
    tuple val(meta), path("*_AMRFinder_Organism.csv"), emit: amrfinder_taxa
    path("versions.yml"),                           emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script = params.ica ? "python ${params.ica_path}/get_taxa_for_amrfinder.py" : "get_taxa_for_amrfinder.py"
    """
    ${script} -t $taxa_file -o ${prefix}_AMRFinder_Organism.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
