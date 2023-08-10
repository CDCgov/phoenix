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
    line=$(head -n1 ${ani_file})
    if [[ "${line}" = "No MASH hit found" ]]; then
        echo "No MASH hit found" > "${meta.id}.fastani.txt"
    else
        def prefix = task.ext.prefix ?: "${meta.id}"
        // terra=true sets paths for bc/wget for terra container paths
        if (params.terra==false) {
            terra = ""
        } else if (params.terra==true) {
            terra = "-t terra"
        } else {
            error "Please set params.terra to either \"true\" or \"false\""
        }
        def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
        """
        db_version=\$(echo ${ani_file} | sed 's/.ani.txt//' | sed 's/${meta.id}_//' )
        # Setup to catch any issues while grabbing date from DB name
        if [[ "\${db_version}" = "" ]]; then
            db_version="REFSEQ_unknown"
        fi
    
        ANI_best_hit_formatter.sh -a $ani_file -n ${prefix} -d \${db_version} $terra
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}
