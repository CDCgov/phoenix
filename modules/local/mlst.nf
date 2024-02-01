process MLST {
    tag "$meta.id"
    label 'process_medium'
    // 2.23.0_07282023 - must edit manually below (line 28)!!!
    container 'quay.io/jvhagey/mlst@sha256:e24cec23ab7300dbb5daf4f4ccd1e4aef9d3befcedae8a14c156cba0bbd952b1'

    input:
    tuple val(meta), path(fasta), val(fairy_outcome), path(taxonomy), path(mlst_db_path)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path("versions.yml")          , emit: versions

    when:
    //if there are scaffolds left after filtering
    "${fairy_outcome[4]}" == "PASSED: More than 0 scaffolds in ${meta.id} after filtering."

    script:
    // helps set correct paths to get database version being used
    if (params.terra==false) { terra = false }
    else if (params.terra==true) { terra = true}
    else { error "Please set params.terra to either \"true\" or \"false\""}
    //define variables
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // mlst is suppose to allow gz and non-gz, but when run in the container (outside of the pipeline) it doesn't work. Also, doesn't work on terra so adding unzip step
    def container = task.container.toString() - "quay.io/jvhagey/mlst@"
    def mlst_version = "2.23.0_01242024"
    def mlst_version_clean = mlst_version.split("_")[0]
    """
    if [[ ${fasta} = *.gz ]]
    then
        unzipped_fasta=\$(basename ${fasta} .gz)
        gunzip --force ${fasta}
    else
        unzipped_fasta=${fasta}
    fi

    if [[ -f ${taxonomy} ]]; then
        genus=\$(head -n7 ${taxonomy} | tail -n1 | cut -d\$'\t' -f2)
        species=\$(head -n8 ${taxonomy} | tail -n1 | cut -d\$'\t' -f2)
    else
        genus="UNKNOWN"
        species="UNKNOWN"
    fi

    mlst --threads $task.cpus \\
        \$unzipped_fasta \\
        > ${prefix}.tsv

    scheme=\$(tail -n1 | cut -d \$'\t' -f2 ${prefix}.tsv)
    if [[ \$scheme == "abaumannii_2" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme abaumannii --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "abaumannii" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme abaumannii_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "ecoli" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme ecoli_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "ecoli_2" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme ecoli --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    # New as of MLST 2.23.0, correctness score update results in some other species outperforming instrinsic ecoli_2 alleles in some cases. Force ecoli to run if ANI taxonomy says so
    elif [[ \${genus,,} == "escherichia" ]]; then
        if [[ \$scheme == "aeromonas" ]]; then
            mv ${prefix}.tsv ${prefix}.OLD-tsv
            mlst --scheme ecoli --threads $task.cpus \$unzipped_fasta > ${prefix}_1.tsv
            mlst --scheme ecoli_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
            cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
            rm ${prefix}_*.tsv
        elif [[ \$scheme == "cfreundii" ]]; then
            mv ${prefix}.tsv ${prefix}.OLD-tsv
            mlst --scheme ecoli --threads $task.cpus \$unzipped_fasta > ${prefix}_1.tsv
            mlst --scheme ecoli_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
            cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
            rm ${prefix}_*.tsv
        elif [[ \$scheme == "senterica" ]]; then
            mv ${prefix}.tsv ${prefix}.OLD-tsv
            mlst --scheme ecoli --threads $task.cpus \$unzipped_fasta > ${prefix}_1.tsv
            mlst --scheme ecoli_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
            cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
            rm ${prefix}_*.tsv
        fi
    else
        :
    fi

    # Add in generic header
    sed -i '1i source_file  Database  ST  locus_1 locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 locus_8 lous_9  locus_10' ${prefix}.tsv

    #handling to get database version being used
    if [[ $terra == false ]]; then
        db_version=\$(cat /mlst-${mlst_version_clean}/db/db_version | date -f - +%Y-%m-%d )
    else
        db_version=\$(cat /opt/conda/envs/phoenix/db/db_version | date -f - +%Y-%m-%d )
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
        mlst_db: \$db_version
        mlst_container: ${container}
    END_VERSIONS
    """
}
