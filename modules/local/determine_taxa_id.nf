process DETERMINE_TAXA_ID {
    tag "${meta.id}"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 25)!!!
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta), path(kraken_weighted), path(formatted_ani), path(k2_bh_summary)
    path(nodes_file)
    path(names_file)

    output:
    tuple val(meta), path('*.tax'), emit: taxonomy
    path("versions.yml"),           emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "bash ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -r needs to be last as in -entry SCAFFOLDS/CDC_SCAFFOLDS k2_bh_summary is not passed so its a blank argument
    def k2_bh_file         = k2_bh_summary ? "-r $k2_bh_summary" : ""
    def k2_weighted_file   = kraken_weighted ? "-k $kraken_weighted" : ""
    def formatted_ani_file = formatted_ani ? "-f $formatted_ani" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}determine_taxID.sh $k2_weighted_file -s ${meta.id} $formatted_ani_file -d $nodes_file -m $names_file $k2_bh_file

    script_version=\$(${ica}determine_taxID.sh -V)

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