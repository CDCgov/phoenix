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
    path(full_outdir)
    val(extended_qc) // true is for -entry CDC_PHOENIX and CDC_SCAFFOLDS

    output:
    tuple val(meta), path('*.scaffolds.fa.gz')    , optional:true, emit: scaffolds // possible that contigs could be created, but not scaffolds
    tuple val(meta), path('*.contigs.fa.gz')      ,                emit: contigs // minimum to complete sucessfully
    tuple val(meta), path('*.assembly.gfa.gz')    , optional:true, emit: gfa
    tuple val(meta), path('*.log')                ,                emit: log
    tuple val(meta), path('*_summaryline.tsv')    , optional:true, emit: line_summary
    tuple val(meta), path('*.synopsis')           , optional:true, emit: synopsis
    path("versions.yml")                          ,                emit: versions
    tuple val(meta), path("*_spades_outcome.csv") ,                emit: spades_outcome

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "bash ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
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
    # Overwrite default that spades was successful
    # Lets downstream process know that spades completed ok - see spades_failure.nf subworkflow
    spades_complete=run_completed
    echo \$spades_complete | tr -d "\\n" > ${prefix}_spades_outcome.csv

    {
        if [[ -z \$(zcat $unpaired_reads) ]]; then
            spades.py \\
                $args \\
                --threads $task.cpus \\
                --memory $maxmem \\
                $input_reads \\
                --phred-offset $phred_offset\\
                -o ./
        else
            spades.py \\
                $args \\
                --threads $task.cpus \\
                --memory $maxmem \\
                $single_reads \\
                $input_reads \\
                --phred-offset $phred_offset\\
                -o ./
        fi

    } || {
    
        if [[ ! -e ${prefix}.contigs.fa.gz ]] || [[ -e ${prefix}.scaffolds.fa.gz ]]; then

            # Set default to be that spades fails and doesn't create scaffolds or contigs
            spades_complete=run_failure
            echo \$spades_complete | tr -d "\\n" > ${prefix}_spades_outcome.csv

            # preemptively create _summary_line.csv and .synopsis file in case spades fails (no contigs or scaffolds created) we can still collect upstream stats. 
            ${ica}pipeline_stats_writer_trimd.sh -a ${fastp_raw_qc} -b ${fastp_total_qc} -c ${reads[0]} -d ${reads[1]} -e ${kraken2_trimd_report} -f ${k2_bh_summary} -g ${krona_trimd}
            ${ica}beforeSpades.sh -k ${k2_bh_summary} -n ${prefix} -d ${full_outdir} ${extended_qc_arg}
        fi
    }

    mv spades.log ${prefix}.spades.log

    #This file will determine if downstream process GENERATE_PIPELINE_STATS_FAILURE and CREATE_SUMMARY_LINE_FAILURE will run (if spades creates contigs, but not scaffolds).
    ${ica}afterSpades.sh

    #get version information
    bspades_version=\$(${ica}beforeSpades.sh -V)
    pipestats_version=\$(${ica}pipeline_stats_writer_trimd.sh -V)
    aspades_version=\$(${ica}afterSpades.sh -V)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
        spades_container: ${container}
        \${bspades_version}
        \${aspades_version}
        \${pipestats_version}
    END_VERSIONS
    """
}
