process CREATE_SAMPLESHEET {
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 23)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    path(directory)
    tuple val(meta), path(assemblies) // LR brings in all assemblies from MEDAKA

    output:
    path("GRiPHin_samplesheet_created.csv"), emit: samplesheet
    path("Assembly_samplesheet.csv"),        emit: assembly_samplesheet
    path("versions.yml"),                    emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/griphin/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def assembly_command = assemblies ? "true" : "false"
    def directory_command = directory ? "true" : "false"
    """

    if [ ${directory_command} = "true" ]; then
        echo "Creating directory samplesheet"
        # creates sample,directory samplesheet for griphin
        ${ica}create_dir_samplesheet.py --directory ${directory}
    elif [ ${assembly_command} = "true" ]; then
        #check if there are fasta files that need to be moved to a folder
        for file in *_medaka_consensus.fasta.gz; do
            if [[ -e "\$file" ]]; then
                mkdir -p assembly_folder
                mv *_medaka_consensus.fasta.gz ./assembly_folder
                break
            fi
        done
        echo "Creating assembly samplesheet"
        # creates sample,assembly samplesheet for griphin
        ${ica}create_assembly_samplesheet.py --assembly assembly_folder
    else
        echo "No valid check type provided, exiting."
        exit 1
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        create_samplesheet.py: \$(${ica}create_samplesheet.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}