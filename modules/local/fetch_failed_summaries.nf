process FETCH_FAILED_SUMMARIES {
    label 'process_single'
    // base_v2.3.0 - MUST manually change below (line 16)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(failed_summaries)
    //path(summaries)

    output:
    path('*_summaryline.tsv'), emit: spades_failure_summary_line
    path("versions.yml")     , emit: versions

    script:
    def container_version = "base_v2.3.0"
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