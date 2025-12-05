process DETERMINE_TOP_MASH_HITS {
    tag "$meta.id"
    label 'process_low'
    // base_v2.3.0 - MUST manually change below (line 31)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(mash_dists), path(assembly_scaffolds)

    output:
    tuple val(meta), path('*_best_MASH_hits.txt'),      emit: top_taxa_list
    tuple val(meta), path('reference_dir'),             emit: reference_dir
    tuple val(meta), path('*_genome_download_log.txt'), emit: log
    path("versions.yml"),                               emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // terra=true sets paths for bc/wget for terra container paths
    def terra = params.terra ? "-t terra" : ""
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "bash ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_name = "${mash_dists}" - ".txt" //get full sample name with REFSEQ_DATE
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    mkdir reference_dir

    ${ica}sort_and_prep_dist.sh -a $assembly_scaffolds -x $mash_dists -o reference_dir $terra

    script_version=\$(${ica}sort_and_prep_dist.sh -V)

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
