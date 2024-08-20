process COLLECT_PROJECT_FILES {
    tag "${meta.id}"
    stageInMode 'copy' // default is symlink. if its not set to copy changes in this script then changes original files.
    label 'process_low'
    container 'quay.io/jvhagey/phoenix:base_v2.1.0'

    input:
    tuple val(meta), path(griphin_excel), path(griphin_tsv), path(phoenix_tsv), path(pipeline_info)

    output:
    tuple val(meta), path("*_old_GRiPHin.tsv"),       emit: griphin_tsv
    tuple val(meta), path("*_old_GRiPHin.xlsx"),      emit: griphin_excel
    tuple val(meta), path("Phoenix_Summary_Old.tsv"), emit: phoenix_tsv
    tuple val(meta), path("software_versions.yml"),   emit: software_versions_file
    path("versions.yml"),                             emit: versions

    script: 
    // define variables
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    mv ./pipeline_info/software_versions.yml .
    #remove the "summary" part of the name so that the output can be picked up correctly in the UPDATE_GRIPHIN process
    mv ${meta.project_id}_GRiPHin_Summary.xlsx ${meta.project_id}_old_GRiPHin.xlsx
    mv ${meta.project_id}_GRiPHin_Summary.tsv ${meta.project_id}_old_GRiPHin.tsv
    mv Phoenix_Summary.tsv Phoenix_Summary_Old.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       phoenix_base_container: ${container}
    END_VERSIONS
    """
}