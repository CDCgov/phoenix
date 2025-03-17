process AMRFINDERPLUS_UPDATE {
    tag "update"
    label 'process_low'
    // 4.0.19-2024-12-18.1
    container 'staphb/ncbi-amrfinderplus@sha256:c257c26454748a797bfa0fbc135f42f1e8d78c68e3f20ba4df804eb234ac9377'

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
