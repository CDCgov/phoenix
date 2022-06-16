process AR_REPORT {
    tag "${meta.id}"
    label 'process_low'

    container 'quay.io/wslh-bioinformatics/ar-report:latest'

    input:
    tuple val(meta), path(ar_report)

    output:
    path "*.ar-report.html"                                      , optional:true,        emit: ar_html
    path "*.ar-report.docx"                                      , optional:true,        emit: ar_docx

    script:
    def report = "--artable ${ar_report}"
    """
    render_report.R \\
    ${report}

    """
}
