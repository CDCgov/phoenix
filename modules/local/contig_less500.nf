process BBMAP_REFORMAT {
    tag "$meta.id"
    label 'process_medium'
    //v39.01
    container 'staphb/bbtools@sha256:161b0e1e198110b7edff8084ae9854d84eb32789d0fd62c7ced302078911c9d7'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*filtered.scaffolds.fa.gz')   , emit: filtered_scaffolds
    tuple val(meta), path('*.log')                       , emit: log
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def raw      = "in=${reads[0]}"
    def trimmed  = "out=${prefix}.filtered.scaffolds.fa.gz"
    def minlength = params.minlength
    def maxmem = task.memory.toGiga()-(task.attempt*9) // keep heap mem low so and rest of mem is for java expansion.
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
    reformat.sh \\
        -Xmx\$maxmem \\
        $raw \\
        $trimmed \\
        threads=$task.cpus \\
        $args \\
        minlength=$minlength \\
        &> ${prefix}.bbmap_filtered.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bbmap: \$(bbversion.sh)
        bbmap_container: ${container}
    END_VERSIONS

    #revert path back to main envs for running on terra
    $terra_exit
    """
}
