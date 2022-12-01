process SCAFFOLDS_SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process-low'
    container 'quay.io/jvhagey/phoenix:base_v1.0.0'

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/fhcrc-microbiome/python-pandas' }"
    
    input:
    path samplesheet
    
    output:
    path '*.valid.csv'       , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """
    scaffolds_samplesheet.py \\
    $samplesheet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
}
