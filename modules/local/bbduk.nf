process BBDUK {
    tag "$meta.id"
    label 'process_medium'
    //v39.13
    container 'staphb/bbtools@sha256:161b0e1e198110b7edff8084ae9854d84eb32789d0fd62c7ced302078911c9d7'

    input:
    tuple val(meta), path(reads)
    path(contaminants)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = meta.single_end ? "in=${reads[0]}" : "in1=${reads[0]} in2=${reads[1]}"
    def trimmed  = meta.single_end ? "out=${prefix}.fastq.gz" : "out1=${prefix}_cleaned_1.fastq.gz out2=${prefix}_cleaned_2.fastq.gz"
    def contaminants_fa = contaminants ? "ref=$contaminants" : ''
    def maxmem = task.memory.toGiga()-(task.attempt*12) // keep heap mem low so and rest of mem is for java expansion.
    def container = task.container.toString() - "staphb/bbtools@"
    //set up for terra
    if (params.terra==false) {
        terra = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra = "PATH=/opt/conda/envs/bbmap/bin:\$PATH"
        terra_exit = """PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/bbmap/bin:||')" """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    """
    #adding path for running bbmap on terra
    $terra

    maxmem=\$(echo \"$maxmem GB\"| sed 's/ GB/g/g')
    bbduk.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $args \\
        $contaminants_fa \\
        &> ${prefix}.bbduk.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
        bbmap_container: ${container}
    END_VERSIONS

    #revert path back to main envs for running on terra
    $terra_exit
    """
}