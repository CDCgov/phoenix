process FORMAT_ANI {
    tag "$meta.id"
    label 'process_single'
    // base_v2.1.0 - MUST manually change below (line 25)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f0304fe170ee359efd2073dcdb4666dddb96ea0b79441b1d2cb1ddc794de4943'

    input:
    tuple val(meta), path(ani_file)

    output:
    tuple val(meta), path('*.fastANI.txt'), emit: ani_best_hit
    path("versions.yml"),                   emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.1.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def script = params.ica ? "${params.ica_path}/ANI_best_hit_formatter.sh" : "ANI_best_hit_formatter.sh"
    def terra = params.terra ? "-t terra" : ""
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
        ${script} -a ${ani_file} -n ${prefix} -d \${db_version} ${terra}
    fi

    script_version=\$(${script} -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}
