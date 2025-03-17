process FETCH_FAILED_SUMMARIES {
    label 'process_single'
    // base_v2.2.0 - MUST manually change below (line 16)!!!
    container 'quay.io/jvhagey/phoenix@sha256:caa2a5660c73d0376d7beb14069436a0e2403bda68904ff140cb789bf4f8753d'

    input:
    path(directory)
    path(failed_summaries)
    path(summaries)

    output:
    path('*_summaryline.tsv'), emit: spades_failure_summary_line
    path("versions.yml")     , emit: versions

    script:
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    #for each summaryline_failure.tsv file check to see if 'SPAdes_Failure' is in the file.
    if [ -f ${directory}/*/*_summaryline_failure.tsv ]; then
        for file in ${directory}/*/*_summaryline_failure.tsv; do 
            if grep -q SPAdes_Failure "\$file"; then
                # if so then add the sample name to the front of the file and move it to the correct place.
                fname=\$(basename \$file _summaryline_failure.tsv)
                cp \$file \${fname}_summaryline.tsv
                mv \$file ${directory}/\${fname}/\${fname}_summaryline.tsv
            fi
        done
    # If the summarylines file doesn't exist then just create an empty file. 
    else
        touch empty_summaryline.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}