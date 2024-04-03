process DETERMINE_TOP_MASH_HITS {
    tag "$meta.id"
    label 'process_low'
    // base_v2.1.0 - MUST manually change below (line 31)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

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
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = "${mash_dists}" - ".txt" //get full sample name with REFSEQ_DATE
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script = params.ica ? "${params.ica_path}/sort_and_prep_dist.sh" : "sort_and_prep_dist.sh"
    def terra = params.terra ? "-t terra" : ""
    """
    mkdir reference_dir

    ${script} -a $assembly_scaffolds -x $mash_dists -o reference_dir $terra

    script_version=\$(${script} -V)

    if [[ ! -f ${sample_name}_best_MASH_hits.txt ]]; then
        echo "No MASH hit found" > ${sample_name}_best_MASH_hits.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Date_of_RefSeq_Pull: \$(date +"%Y-%m-%d")
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}
