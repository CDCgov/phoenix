process FORMAT_ANI {
    tag "$meta.id"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 25)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    tuple val(meta), path(ani_file)

    output:
    tuple val(meta), path('*.fastANI.txt'), optional:true, emit: ani_best_hit
    tuple val(meta), path('*.to_check_fastANI.txt'),       emit: ani_best_hit_to_check
    path("versions.yml"),                                  emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // terra=true sets paths for bc/wget for terra container paths
    def terra = params.terra ? "-t terra" : ""
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
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

        # since we need to check any files that have Escherichia or Shigella in them we will rename files
        if grep -qE "Escherichia|Shigella" "${prefix}_\${db_version}.fastANI_initial.txt"; then
            mv ${prefix}_\${db_version}.fastANI_initial.txt ${prefix}_\${db_version}.to_check_fastANI.txt
        else
            cp ${prefix}_\${db_version}.fastANI_initial.txt ${prefix}_\${db_version}.to_check_fastANI.txt
            mv ${prefix}_\${db_version}.fastANI_initial.txt ${prefix}_\${db_version}.fastANI.txt
        fi
    fi

    script_version=\$(${ica}ANI_best_hit_formatter.sh -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
        \${script_version}
    END_VERSIONS
    """
}
