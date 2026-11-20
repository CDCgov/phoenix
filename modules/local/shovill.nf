process SHOVILL {
    tag "$meta.id"
    label 'process_high_memory'
    // shovill v1.4.2 with spades v3.15.5
    container 'staphb/shovill@sha256:68d58389830ccce6057b038efcd6a70219253f4ffee945c47a98998a3da774ca'
    
    input:
    tuple val(meta), path(reads), path(unpaired_reads), path(k2_bh_summary),
    path(fastp_raw_qc),
    path(fastp_total_qc), 
    path(kraken2_trimd_report),
    path(krona_trimd), 
    path(full_outdir),
    path(outcome_to_edit)
    val(extended_qc) // true is for -entry CDC_PHOENIX and CDC_SCAFFOLDS

    output:
    tuple val(meta), path('*.scaffolds.fa.gz'),        optional:true, emit: scaffolds
    tuple val(meta), path('*.contigs.fa.gz'),          optional:true, emit: contigs
    tuple val(meta), path('*_trimstats_summary.txt'),  optional:true, emit: fairy_outcome
    tuple val(meta), path('*.assembly.gfa.gz'),        optional:true, emit: gfa
    tuple val(meta), path('*.log'),                    emit: log
    path('versions.yml'),                              emit: versions
    // Keeping SPAdes outcome naming for compatibility
    tuple val(meta), path('*_spades_outcome.csv'),     emit: spades_outcome

    script:
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
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
    #adding shovill path for running shovill on terra
    $terra

    # Default status: assembly failure
    printf 'run_failure' > ${prefix}_spades_outcome.csv

    # Run Shovill using its SPAdes backend
    if shovill \\
        --R1 ${reads[0]} \\
        --R2 ${reads[1]} \\
        --outdir shovill_output \\
        --assembler spades \\
        --cpus ${task.cpus} \\
        --ram ${maxmem} \\
        --minlen 500 \\
        --namefmt ${prefix}_%05d \\
        ${args} \\
        > ${prefix}.shovill.log 2>&1
    then
        printf 'run_completed' > ${prefix}_spades_outcome.csv
    else
        echo "Shovill failed for ${prefix}; preserving failure status"
    fi

    # Normalize Shovill outputs for PHoeNIx
    ${ica}afterShovill.py \\
        --prefix ${prefix} \\
        --outdir shovill_output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$(shovill --version 2>&1 | sed 's/^.*shovill //; s/ .*\$//')
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
        shovill_container: ${container}
        \$(${ica}afterShovill.py -V)
    END_VERSIONS

    #revert path back to main envs for running on terra
    $terra_exit
    """
}
