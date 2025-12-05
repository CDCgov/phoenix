process RENAME_FASTA_HEADERS {
    tag "$meta.id"
    label 'process_low'
    // base_v2.3.0 - MUST manually change below (line 21)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    tuple val(meta), path(assembled_scaffolds)

    output:
    tuple val(meta), path('*.renamed.scaffolds.fa.gz'), emit: renamed_scaffolds
    path "versions.yml"                               , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def prefix = task.ext.prefix ?: "${meta.id}"
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    """
    gunzip --force ${assembled_scaffolds}
    unzipped=\$(basename ${assembled_scaffolds} .gz) #adding this in to allow alternative file names with -entry SCAFFOLDS --scaffolds_ext

    ${ica}rename_fasta_headers.py --input \$unzipped --output ${prefix}.renamed.scaffolds.fa --name ${prefix}

    gzip --force ${prefix}.renamed.scaffolds.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        rename_fasta_headers.py: \$(${ica}rename_fasta_headers.py --version )
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}