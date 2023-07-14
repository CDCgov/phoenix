process DETERMINE_TAXA_ID {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.0.0'

    input:
    tuple val(meta), path(kraken_weighted), path(formatted_ani_file), path(k2_bh_summary)
    path(taxa_file)

    output:
    tuple val(meta), path('*.tax'), emit: taxonomy
    path("versions.yml")          , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -r needs to be last as in -entry SCAFFOLDS/CDC_SCAFFOLDS k2_bh_summary is not passed so its a blank argument
    def k2_bh_file = k2_bh_summary ? "-r $k2_bh_summary" : ""
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    determine_taxID.sh -k $kraken_weighted -s $meta.id -f $formatted_ani_file -d $taxa_file $k2_bh_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NCBI Taxonomy Reference File: $taxa_file
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}