process MLST {
    tag "$meta.id"
    label 'process_medium'
    // 2.25.0_12102025 - must edit manually below (line 28)!!!
    container 'quay.io/jvhagey/mlst@sha256:75da09f926c18c69c656ff22c5bbee6f037307aa4a41d11b14ff0ccaa2d6c628'

    input:
    tuple val(meta), path(fasta), path(taxonomy)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path("versions.yml")          , emit: versions

    script:
    // helps set correct paths to get database version being used
    def terra = params.terra ? "true" : "false"
    //define variables
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // mlst is suppose to allow gz and non-gz, but when run in the container (outside of the pipeline) it doesn't work. Also, doesn't work on terra so adding unzip step
    def container = task.container.toString() - "quay.io/jvhagey/mlst@"
    def mlst_version = "2.25.0_12102025"
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
    elif [[ \${genus,,} == "mycobacterium" ]]; then
        if [[ \$scheme == "mabscessus" ]]; then
            mv ${prefix}.tsv ${prefix}_1.tsv
            mlst --scheme mycobacteria --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
            cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
            rm ${prefix}_*.tsv
        elif [[ \$scheme == "mycobacteria" && \${species,,} == "abscessus" ]]; then
            mv ${prefix}.tsv ${prefix}_1.tsv
            mlst --scheme mabscessus --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
            cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
            rm ${prefix}_*.tsv
        fi

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
    elif [[ \${genus,,} == "enterobacter" ]]; then
        if [[ \$scheme == "cronobacter" ]]; then
            mv ${prefix}.tsv ${prefix}.OLD-tsv
            mlst --scheme ecloacae --threads $task.cpus \$unzipped_fasta > ${prefix}.tsv
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
