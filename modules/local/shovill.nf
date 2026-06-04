process SHOVILL {
    tag "$meta.id"
    label 'process_high_memory'
    // shovill v1.1.0 with spades v3.15.5
    container 'staphb/shovill@sha256:68f2380df89e5ac6e1bb1e14a03a6589f2dfd0cfdb5e47ca1e3f5810e8b7ce4c'

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
    tuple val(meta), path('*.scaffolds.fa.gz')        , optional:true, emit: scaffolds
    tuple val(meta), path('*.contigs.fa.gz')          , optional:true, emit: contigs
    tuple val(meta), path('*_trimstats_summary.txt')  , optional:true, emit: fairy_outcome
    tuple val(meta), path('*.assembly.gfa.gz')        , optional:true, emit: gfa
    tuple val(meta), path('*.log')                    ,                emit: log
    path("versions.yml")                              ,                emit: versions
    tuple val(meta), path("*_spades_outcome.csv")     ,                emit: spades_outcome

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "bash ${params.bin_dir}" : ""
    // define variables
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def container = task.container.toString() - "staphb/shovill@"
    //set up for terra
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = "PATH=/opt/conda/envs/shovill/bin:\$PATH"
        terra_exit = """PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/shovill/bin:||')" """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
    #adding shovill path for running on terra
    $terra

    # Set default to be that shovill fails and doesn't create scaffolds or contigs
    spades_complete=run_failure
    echo \$spades_complete | tr -d "\\n" > ${prefix}_spades_outcome.csv

    {
        shovill \\
            --outdir shovill_output \\
            --R1 ${reads[0]} \\
            --R2 ${reads[1]} \\
            --cpus $task.cpus \\
            --ram $maxmem \\
            --minlen ${params.minlength ?: 500} \\
            --mincov ${params.mincov ?: 2} \\
            --assembler spades \\
            --force \\
            $args

        # Copy shovill outputs to match what afterSpades.sh expects (SPAdes-style names in working dir)
        if [ -f shovill_output/contigs.fa ]; then
            # Shovill's contigs.fa = filtered scaffolds, place as both scaffolds.fasta and contigs.fasta
            cp shovill_output/contigs.fa scaffolds.fasta
            cp shovill_output/contigs.fa contigs.fasta
        fi

        if [ -f shovill_output/contigs.gfa ]; then
            cp shovill_output/contigs.gfa assembly_graph_with_scaffolds.gfa
        fi

        # Overwrite default - shovill was successful
        spades_complete=run_completed
        echo \$spades_complete | tr -d "\\n" > ${prefix}_spades_outcome.csv

    } || {
        # Shovill failed - outcome file already set to failure above
        echo "Shovill failed, keeping failure status and summary files"
    }

    # Create log file matching *.spades.log pattern that afterSpades.sh expects
    if [ -f shovill_output/shovill.log ]; then
        cp shovill_output/shovill.log ${prefix}.spades.log
    else
        echo "Shovill did not produce a log file" > ${prefix}.spades.log
    fi

    #This file will determine if downstream process GENERATE_PIPELINE_STATS_FAILURE and CREATE_SUMMARY_LINE_FAILURE will run
    ${ica}afterSpades.sh

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$(shovill --version 2>&1 | sed 's/^.*shovill //')
        shovill_container: ${container}
        \$(${ica}afterSpades.sh -V)
    END_VERSIONS

    #removing spades/shovill path for running on terra
    $terra_exit
    """
}
