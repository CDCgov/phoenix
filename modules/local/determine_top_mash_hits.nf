process DETERMINE_TOP_MASH_HITS {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(mash_dists), path(assembly_scaffolds), val(fairy_outcome)

    output:
    tuple val(meta), path('*_best_MASH_hits.txt'), emit: top_taxa_list
    tuple val(meta), path('reference_dir'),        emit: reference_dir
    path("versions.yml"),                          emit: versions

    when:
    //if there are scaffolds left after filtering
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) { terra = ""} 
    else if (params.terra==true) { terra = "-t terra" }
    else { error "Please set params.terra to either \"true\" or \"false\"" }
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = "${mash_dists}" - ".txt" //get full sample name with REFSEQ_DATE
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    mkdir reference_dir

    ${ica}sort_and_prep_dist.sh -a $assembly_scaffolds -x $mash_dists -o reference_dir $terra

    if [[ ! -f ${sample_name}_best_MASH_hits.txt ]]; then
        echo "No MASH hit found" > ${sample_name}_best_MASH_hits.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Date of RefSeq Pull: \$(date +"%d-%m-%y")
        phoenix_base_container: ${container}
        \$(${ica}sort_and_prep_dist.sh -V)
    END_VERSIONS
    """
}
