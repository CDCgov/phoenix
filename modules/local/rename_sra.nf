process RENAME_SRA_FASTA {
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rename:1.601--hdfd78af_1' :
        'quay.io/biocontainers/rename:1.601--hdfd78af_1' }"

    input:
    path fastqs
    path version_file

    output:
    path "versions.yml"                                   , emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """ 
    cd $fastqs
    rename 's/_1.fastq.gz/_R1_001.fastq.gz/' *
    rename 's/_2.fastq.gz/_R2_001.fastq.gz/' *

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqtk: \$(echo \$(rename 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    
    """
}