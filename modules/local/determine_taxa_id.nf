process DETERMINE_TAXA_ID {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(kraken_weighted), path(formatted_ani_file), path(k2_bh_summary)
    path(nodes_file)
    path(names_file)

    output:
    tuple val(meta), path('*.tax'), emit: taxonomy
    path("versions.yml")          , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -r needs to be last as in -entry SCAFFOLDS/CDC_SCAFFOLDS k2_bh_summary is not passed so its a blank argument
    def k2_bh_file = k2_bh_summary ? "-r $k2_bh_summary" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """

    ${ica}determine_taxID.sh -k $kraken_weighted -s $meta.id -f $formatted_ani_file -d $nodes_file -m $names_file $k2_bh_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NCBI Taxonomy Nodes Reference File: $nodes_file
        NCBI Taxonomy Names Reference File: $names_file
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}