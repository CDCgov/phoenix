process PLSDB {
    tag "$meta.id"
    label 'process_high'
    container 'quay.io/aharring83/plsdb@sha256:6674c718325ba52174c6309b480d449ced2abf80b95e110e7edabada984ab217' 

    input:
    tuple val(meta), path(fasta)
    //path plsdb
    path plsdb_dir
    path conf
    

    output:
    tuple val(meta), path("${meta.id}_plasmidID.tsv"), emit: plasmidID
    tuple val(meta), path("${meta.id}_contigs/*"), emit: contigs
    path ("versions.yml"), emit: versions
    
    script:

    def container_version = "v0.1.0"//Update this version number when the script is updated
    def container = task.container.toString() - "quay.io/aharring83/plsdb:"
    def plsdb_version = "20240521.v2" //Update this version number when the database is updated

    """
    plasmid.py -i $fasta -db $plsdb_dir/plsdb_2024_05_31_v2 -c $conf -o ${meta.id}_plasmidID.tsv
    mv contigs ${meta.id}_contigs
    #if [ -f data.json ]; then
    #cat ${meta.id}_plasmid.tsv  
    #else
    #touch ${meta.id}_plasmid.tsv
    #fi
    
    # NEED TO UPDATE THE PLSDB VERSION NUMBER WHEN THE DATABASE IS UPDATED/MAKE NEW VERSION
    # handle the plsdb container with versioning later, but for now just hardcoding the version number in the versions.yml file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        plasmid.py: \$(plasmid.py --version)
        #plsdb_creation_date: 20240521.v2
        plsdb_container: ${container}
        plsdb_container_version: ${container_version}
        plsdb_version: ${plsdb_version}
        
    END_VERSIONS

    """

}