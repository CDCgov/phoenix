process FAIRY {
    tag "$meta.id"
    label 'process_low'
    label 'error_ignore'
    container 'quay.io/jvhagey/phoenix:base_v1.1.0'

    input:
    tuple val(meta), path(reads)

    output:
    
    tuple val(meta), path('*_results.txt')       , emit: results
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fname = "${reads[0]}".minus(".fastq.gz")
    def fnameB = "${reads[1]}".minus(".fastq.gz")

    """
    gzip -t ${reads[0]} 2>> ${fname}.txt
    gzip -t ${reads[1]} 2>> ${fnameB}.txt
    
    if grep -Fxq "error" ${fname}.txt || grep -Fxq "error" ${fnameB}.txt; then
        echo "FAIL" > ${prefix}_results.txt
    else
        echo "PASS" > ${prefix}_results.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}