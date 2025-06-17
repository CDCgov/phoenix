process COLLECT_PROJECT_FILES {
    tag "${meta.project_id}"
    stageInMode 'copy' // default is symlink. if its not set to copy changes in this script then changes original files.
    label 'process_low'
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    tuple val(meta), path(griphin_excel), path(griphin_tsv), path(phoenix_tsv), path(pipeline_info)
    val(full_project_dir)

    output:
    tuple val(meta), path("*_old_GRiPHin.tsv"),       emit: griphin_tsv
    tuple val(meta), path("*_old_GRiPHin.xlsx"),      emit: griphin_excel
    tuple val(meta), path("Phoenix_Summary_Old.tsv"), emit: phoenix_tsv
    tuple val(meta), path("software_versions.yml"),   emit: software_versions_file
    path("versions.yml"),                             emit: versions

    script: 
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def parentPath = meta.project_id.tokenize('/').last()
    def project_dir = full_project_dir ? "${parentPath}" : "${meta.project_id}"
    """
    if [ ! -e "software_versions.yml" ]; then
        mv ./pipeline_info/software_versions.yml .
    fi
    #remove the "summary" part of the name so that the output can be picked up correctly in the UPDATE_GRIPHIN process
    mv ${project_dir}_GRiPHin_Summary.xlsx ${project_dir}_old_GRiPHin.xlsx
    mv ${project_dir}_GRiPHin_Summary.tsv ${project_dir}_old_GRiPHin.tsv
    mv Phoenix_Summary.tsv Phoenix_Summary_Old.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}