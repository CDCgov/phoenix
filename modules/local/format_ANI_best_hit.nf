process FORMAT_ANI {
    tag "$meta.id"
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(ani_file)

    output:
    tuple val(meta), path('*.fastANI.txt'), emit: ani_best_hit
    path("versions.yml"),                   emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
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
    ANI_best_hit_formatter.sh -a $ani_file -n ${prefix} $terra

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}