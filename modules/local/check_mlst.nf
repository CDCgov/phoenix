process CHECK_MLST {
    tag "$meta.id"
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 23)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    tuple val(meta), path(mlst_file), path(taxonomy_file), path(local_dbases)

    output:
    tuple val(meta), path("*_combined.tsv"), emit: checked_MLSTs
    tuple val(meta), path("*_status.txt"),   emit: status
    path("versions.yml")                 ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    ${ica}fix_MLST2.py --input $mlst_file --taxonomy $taxonomy_file --mlst_database ${local_dbases}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        fix_MLST2.py: \$(${ica}fix_MLST2.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}