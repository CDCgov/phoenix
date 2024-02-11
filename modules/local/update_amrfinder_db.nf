process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_low'
    // 3.11.11-2023-04-17.1
    container 'staphb/ncbi-amrfinderplus@sha256:194eec0c758f92c3c8a8884b9f1ddbfb7626977459e8938a6ece98aceb8e3bbd'

    output:
    path "amrfinderdb.tar.gz", emit: db
    path "versions.yml"      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def container = task.container.toString() - "quay.io/jvhagey/ncbi-amrfinderplus@"
    """
    mkdir amrfinderdb
    amrfinder_update -d amrfinderdb
    version=\$(head amrfinderdb/latest/version.txt)
    tar czvf amrfinderdb.tar.gz -C ./amrfinderdb/\$version ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus_db_version: \$(head amrfinderdb/latest/version.txt)
        amrfinderplus_container: ${container}
    END_VERSIONS
    """
}
