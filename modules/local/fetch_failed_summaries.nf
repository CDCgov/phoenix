process FETCH_FAILED_SUMMARIES {
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 16)!!!
    container 'quay.io/jvhagey/phoenix@sha256:f7cb3aa4e3324cab43d8635be17da8ae15f62e39d380acda844d1c9deef69c60'

    input:
    tuple val(meta), path(failed_summaries)
    //path(summaries)

    output:
    path('*_summaryline.tsv'), emit: spades_failure_summary_line
    path("versions.yml")     , emit: versions

    script:
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    #summaryline_failure.tsv rename and publish
    mv ${meta.id}_summaryline_failure.tsv ${meta.id}_summaryline.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}