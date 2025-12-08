process GET_TAXA_FOR_AMRFINDER {
    tag "$meta.id"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 21)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f7cb3aa4e3324cab43d8635be17da8ae15f62e39d380acda844d1c9deef69c60'

    input:
    tuple val(meta), path(taxa_file)

    output:
    tuple val(meta), path("*_AMRFinder_Organism.csv"),               emit: amrfinder_taxa
    tuple val(meta), path("*_ABRITAMR_Organism.csv"), optional:true, emit: abritamr_taxa
    path("versions.yml"),                                            emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.2.0"
    def abritamr_val = (params.mode_upper == "CLIA" || params.run_abritamr == true) ? "--abritamr_taxa" : "" //add if you are running the clia version to get the taxa for abritamr
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}get_taxa_for_amrfinder.py -t $taxa_file -o ${prefix} ${abritamr_val}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
