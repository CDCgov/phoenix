process CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_low'
    tag "${meta.project_id}"

    // Requires `pyyaml` which does not have a dedicated container but is in the MultiQC container
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(versions)

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    // 📦 This replaces the logic in modules.config. Must use a seperate val var because of runtime restrictions
    //publishDir "${project_dir_val}/centar_pipeline_info", mode: 'copy', pattern: '*_versions.yml'

    script:
    def args = task.ext.args ?: ''
    template 'dumpsoftwareversions.py'
}