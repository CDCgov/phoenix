process GET_TAXA_FOR_AMRFINDER {
    tag "$meta.id"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 21)!!!
    container 'quay.io/jvhagey/phoenix@sha256:caa2a5660c73d0376d7beb14069436a0e2403bda68904ff140cb789bf4f8753d'

    input:
    tuple val(meta), path(taxa_file)
    val(clia_entry)

    output:
    tuple val(meta), path("*_AMRFinder_Organism.csv"),               emit: amrfinder_taxa
    tuple val(meta), path("*_ABRITAMR_Organism.csv"), optional:true, emit: abritamr_taxa
    path("versions.yml"),                                            emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.2.0"
    def clia = clia_entry ? "--clia" : "" //add if you are running the clia version to get the taxa for abritamr
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}get_taxa_for_amrfinder.py -t $taxa_file -o ${prefix} ${clia}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
