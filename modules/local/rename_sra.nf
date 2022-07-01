process RENAME_SRA_FASTA {
    label 'process_low'

    input:
    path versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    """
    rename_sra.py
    
    """
}