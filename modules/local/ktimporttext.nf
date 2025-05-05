process KRONA_KTIMPORTTEXT {
    tag "$meta.id"
    label 'process_single'
    // 2.8.1--pl5321hdfd78af_1
    container 'quay.io/biocontainers/krona@sha256:8917b9840b369d102ee759a037cc8577295875952013aaa18897c00569c9fe47'

    input:
    tuple val(meta), path(krona)
    val(type) //weighted, trimmmed or assembled

    output:
    tuple val(meta), path('*.html'), emit: html
    path("versions.yml")           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container = task.container.toString() - "quay.io/biocontainers/krona@"
    """
    ktImportText  \\
        $args \\
        -o ${prefix}_${type}.html \\
        $krona

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        krona: \$( echo \$(ktImportText 2>&1) | sed 's/^.*KronaTools //g; s/- ktImportText.*\$//g')
        krona_container: ${container}
    END_VERSIONS
    """
}
