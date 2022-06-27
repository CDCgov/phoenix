process SRATOOLS_FASTERQDUMP {
    tag "$sra"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::sra-tools=2.11.0 conda-forge::pigz=2.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:6a9ff0e76ec016c3d0d27e0c0d362339f2d787e6-0' :
        'quay.io/biocontainers/mulled-v2-5f89fe0cd045cb1d615630b9261a1d17943a9b6a:6a9ff0e76ec016c3d0d27e0c0d362339f2d787e6-0' }"

    input:
    path(sra)

    output:
    path(output)                 , emit: reads
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    
    """
    fasterq-dump \\
        $args \\
        --threads $task.cpus \\
        ${sra.name} 

    pigz \\
        $args2 \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}