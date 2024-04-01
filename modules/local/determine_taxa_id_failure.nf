process DETERMINE_TAXA_ID_FAILURE {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 26)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(k2_bh_summary), val(spades_outcome)
    path(nodes_file)
    path(names_file)

    output:
    tuple val(meta), path('*.tax'), emit: taxonomy
    path("versions.yml"),           emit: versions

    when:
    "${spades_outcome[0]}" == "run_failure" || "${spades_outcome[1]}" == "no_scaffolds" || "${spades_outcome[2]}" == "no_contigs"

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script = params.ica ? "bash ${params.ica_path}/determine_taxID.sh" : "determine_taxID.sh"
    """
    ${script} -r $k2_bh_summary -s $meta.id -d $nodes_file -m $names_file

    script_version=\$(${script} -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NCBI_Taxonomy_Nodes_Reference_File: $nodes_file
        NCBI_Taxonomy_Names_Reference_File: $names_file
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}