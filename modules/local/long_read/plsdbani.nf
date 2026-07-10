process PLSDBANI {
    tag "$meta.id"
    label 'process_high'
    //errorStrategy 'ignore'
    //container 'aharring83/rscript'
    container "quay.io/aharring83/rscript@sha256:8b97e21c343a9abfab74e91856ece45c0181a350c586206ac9fc7807949d904f"

    input:
    tuple val(meta), path(plasmidID), path(fasta)
    path(viz)
    path(plsdbfasta)


    output:
    tuple val(meta), path("${meta.id}_plasmidANI.csv"), emit: plasmidANI
    tuple val(meta), path("accession"), emit: plasmidaccession
    tuple val(meta), path("contig"), emit: contigs
    tuple val(meta), path("Viz"), emit: visual
    path ("versions.yml"), emit: versions
    
    script:

    def container_version = "v0.1.0" //Update this version number when the script is updated
    def container = task.container.toString() - "quay.io/aharring83/rscript:"

    """
    ani.py -i $plasmidID -a $fasta 
    
    mkdir Viz
    mv *.visual Viz
    #ani.sh
    awk 'BEGIN {OFS=","; print "query,reference, ANI%,total,aligned"} {print \$1, \$2, \$3, \$4,\$5}' *.txt > ${meta.id}_plasmidANI.csv
    viz.py $plasmidID contig/ accession/ Viz $viz
     
    #mkdir ${meta.id}
    #mv contig ${meta.id}/
    #mv accession ${meta.id}/
    #mv Viz ${meta.id}/
    #find . -maxdepth 1 -type f ! -name 'versions.yml' -exec mv -t ${meta.id} {} +   #move all files except versions.yml to the ${meta.id} directory

    # NEED TO UPDATE THE PLSDB VERSION NUMBER WHEN THE DATABASE IS UPDATED/MAKE NEW VERSION
    # handle the plsdb container with versioning later, but for now just hardcoding the version number in the versions.yml file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        ani.py: \$(ani.py --version)
        viz.py: \$(viz.py --version)
        plsdb_creation_date: 20240521.v2
        rscript_container: ${container}
        
    END_VERSIONS

    """
}