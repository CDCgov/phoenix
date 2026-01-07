process COLLECT_PROJECT_FILES {
    tag "${meta.project_id.toString().split('/')[-1].replace(']', '')}"
    //tag "${meta.project_id}"
    stageInMode 'copy' // default is symlink. if its not set to copy changes in this script then changes original files.
    label 'process_low'
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    tuple val(meta), path(griphin_excel), path(griphin_tsv), path(phoenix_tsv), path(pipeline_info)
    val(full_project_dir)

    output:
    tuple val(meta), path("*_old_GRiPHin.tsv"),       emit: griphin_tsv
    tuple val(meta), path("*_old_GRiPHin.xlsx"),      emit: griphin_excel
    tuple val(meta), path("Phoenix_Summary_Old.tsv"), emit: phoenix_tsv
    tuple val(meta), path("*_software_versions.yml"), emit: software_versions_file
    path("versions.yml"),                             emit: versions

    script: 
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def parentPath = meta.project_id.tokenize('/').last()
    def project_dir = full_project_dir ? "${parentPath}" : "${meta.project_id}"
    """

    if [ ! -e "software_versions.yml" ]; then
        mv ./pipeline_info/software_versions.yml ${project_dir}_software_versions.yml
    else
        mv software_versions.yml ${project_dir}_software_versions.yml
    fi
    #remove the "summary" part of the name so that the output can be picked up correctly in the UPDATE_GRIPHIN process
    mv *_GRiPHin_Summary.xlsx ${project_dir}_old_GRiPHin.xlsx
    mv *_GRiPHin_Summary.tsv ${project_dir}_old_GRiPHin.tsv
    mv Phoenix_Summary.tsv Phoenix_Summary_Old.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}