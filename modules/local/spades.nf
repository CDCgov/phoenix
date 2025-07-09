process SPADES {
    tag "$meta.id"
    label 'process_high_memory'
    // v4.0.0
    container 'staphb/spades@sha256:5df39e8404df2678ccc6c6ed9d7aa0e59b79dfa798aef7fd4fc06cc86ba0b4c0'

    input:
    tuple val(meta), path(reads), path(unpaired_reads), path(k2_bh_summary), \
    path(fastp_raw_qc), \
    path(fastp_total_qc), \
    path(kraken2_trimd_report), \
    path(krona_trimd), \
    path(full_outdir), \
    path(outcome_to_edit)
    val(extended_qc) // true is for -entry CDC_PHOENIX and CDC_SCAFFOLDS

    output:
    tuple val(meta), path('*.scaffolds.fa.gz')        , optional:true, emit: scaffolds // possible that contigs could be created, but not scaffolds
    tuple val(meta), path('*.contigs.fa.gz')          , optional:true, emit: contigs // minimum to complete sucessfully - changed to optional:true to allow for failure handling
    tuple val(meta), path('*_trimstats_summary.txt')  , optional:true, emit: fairy_outcome
    tuple val(meta), path('*.assembly.gfa.gz')        , optional:true, emit: gfa
    tuple val(meta), path('*.log')                    ,                emit: log
    //tuple val(meta), path('*_summaryline_failure.tsv'), optional:true, emit: line_summary_failure
    path("versions.yml")                              ,                emit: versions
    tuple val(meta), path("*_spades_outcome.csv")     ,                emit: spades_outcome

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "bash ${params.bin_dir}" : ""
    // define variables
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga() // allow 4 less GB to provide enough space
    def input_reads = "-1 ${reads[0]} -2 ${reads[1]}"
    def single_reads = "-s $unpaired_reads"
    def phred_offset = params.phred
    def extended_qc_arg = extended_qc ? "-c" : ""
    def container = task.container.toString() - "staphb/spades@"
    """
    # Set default to be that spades fails and doesn't create scaffolds or contigs
    spades_complete=run_failure
    echo \$spades_complete | tr -d "\\n" > ${prefix}_spades_outcome.csv

    {
        if [[ -z \$(zcat $unpaired_reads) ]]; then
            echo "I am here with paired reads"
            spades.py \\
                $args \\
                --threads $task.cpus \\
                --memory $maxmem \\
                $input_reads \\
                --phred-offset $phred_offset\\
                -o ./
        else
            echo "I am here with unpaired reads"
            spades.py \\
                $args \\
                --threads $task.cpus \\
                --memory $maxmem \\
                $single_reads \\
                $input_reads \\
                --phred-offset $phred_offset\\
                -o ./
        fi

        # Overwrite default that spades was successful
        # Lets downstream process know that spades completed ok - see spades_failure.nf subworkflow
        spades_complete=run_completed
        echo \$spades_complete | tr -d "\\n" > ${prefix}_spades_outcome.csv

    } || {
    
        # SPAdes failed - outcome file already set to failure above
        echo "SPAdes failed, keeping failure status and summary files"
    }

    mv spades.log ${prefix}.spades.log

    #This file will determine if downstream process GENERATE_PIPELINE_STATS_FAILURE and CREATE_SUMMARY_LINE_FAILURE will run (if spades creates contigs, but not scaffolds).
    # also, will rename the fairy file to publish
    ${ica}afterSpades.sh

    #get version information
    aspades_version=\$(${ica}afterSpades.sh -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
        spades_container: ${container}
        \${aspades_version}
    END_VERSIONS
    """
}
