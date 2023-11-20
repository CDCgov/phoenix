process FORMAT_ANI {
    tag "$meta.id"
    label 'process_single'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(ani_file)

    output:
    tuple val(meta), path('*.fastANI.txt'), emit: ani_best_hit
    path("versions.yml"),                   emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // terra=true sets paths for bc/wget for terra container paths
    if (params.terra==false) { terra = ""} 
    else if (params.terra==true) { terra = "-t terra" }
    else { error "Please set params.terra to either \"true\" or \"false\"" }
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${workflow.launchDir}/bin/" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    line=\$(head -n1 ${ani_file})
    if [[ "\${line}" == "Mash/FastANI Error:"* ]]; then
        echo "Mash/FastANI Error: No MASH hit found." > "${prefix}.fastANI.txt"
    else
        db_version=\$(echo ${ani_file} | sed 's/.ani.txt//' | sed 's/${prefix}_//' )
        # Setup to catch any issues while grabbing date from DB name
        if [[ "\${db_version}" = "" ]]; then
            db_version="REFSEQ_unknown"
        fi
        # script also checks that match is 80 or > otherwise an error is thrown
        ${ica}ANI_best_hit_formatter.sh -a ${ani_file} -n ${prefix} -d \${db_version} ${terra}
    fi

    script_version=\$(${ica}ANI_best_hit_formatter.sh -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}
