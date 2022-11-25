process SCAFFOLDS_SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process-low'
    
    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/fhcrc-microbiome/python-pandas' }"
}
