process SPADES {
    tag "$meta.id"
    label 'process_high'
    container 'staphb/spades:3.15.4'
    afterScript "sh ${baseDir}/bin/afterSpades.sh" // Handles file zipping, renaming with prefix and checks that files were created
    // Create a summaryline file that will be deleted later if spades is successful if not this line shows up in the final Phoenix_output_summary file
    beforeScript "echo -e '${meta.id}\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tFAIL\tSPAdes_Failure' | tr -d '\n' > ${meta.id}_summaryline_failure.tsv ; cp ${meta.id}_summaryline_failure.tsv ${params.outdir}/${meta.id}"

    input:
    tuple val(meta), path(reads), path(unpaired_reads)

    output:
    tuple val(meta), path('*.scaffolds.fa.gz')    , optional:true, emit: scaffolds
    tuple val(meta), path('*.contigs.fa.gz')      ,                emit: contigs
    tuple val(meta), path('*.assembly.gfa.gz')    , optional:true, emit: gfa
    tuple val(meta), path('*.log')                ,                emit: log
    tuple val(meta), path('*_summaryline.tsv')    , optional:true, emit: line_summary
    path  "versions.yml"                          ,                emit: versions
    tuple val(meta), path("*_spades_outcome.csv") ,                emit: spades_outcome

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def input_reads = "-1 ${reads[0]} -2 ${reads[1]}"
    def single_reads = "-s $unpaired_reads"
    def phred_offset = params.phred
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS

    spades_complete=run_failure,no_scaffolds,no_contigs
    echo \$spades_complete | tr -d "\\n" > ${prefix}_spades_outcome.csv

    spades.py \\
        $args \\
        --threads $task.cpus \\
        --memory $maxmem \\
        $single_reads \\
        $input_reads \\
        --phred-offset $phred_offset\\
        -o ./

    mv spades.log ${prefix}.spades.log
    spades_complete=run_completed
    echo \$spades_complete | tr -d "\\n" > ${prefix}_spades_outcome.csv
    rm ${params.outdir}/${meta.id}/${meta.id}_summaryline_failure.tsv
    """
}