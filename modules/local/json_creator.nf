process JSON_CREATOR {
    tag "$file"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::pandas=1.1.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.1.5' :
        'quay.io/fhcrc-microbiome/python-pandas' }"

    input:
    path file

    output:
    path '*.json'       , emit: json
    path "versions.yml", emit: versions

    script: 
    """
    to_json.py \\
    $file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
