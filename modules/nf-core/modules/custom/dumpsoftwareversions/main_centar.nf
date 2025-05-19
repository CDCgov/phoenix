process CENTAR_CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_low'
    tag "${project_id_dir}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    path versions
    path project_path

//    output:
//    path 'CENTAR_software_versions.yml',        emit: yml
//    path "software_versions_mqc.yml", emit: mqc_yml
//    path "versions.yml",              emit: versions

    script:
    """
    dumpsoftwareversions_mod_with_out.py \
        --versions ${versions} \
        --outdir ${project_path} \
        --workflow_name ${workflow.manifest.name} \
        --workflow_version ${workflow.manifest.version} \
        --nextflow_version ${workflow.nextflow.version} \
        --source "CENTAR"
    """
}