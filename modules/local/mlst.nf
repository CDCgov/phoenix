process MLST {
    tag "$meta.id"
    label 'process_medium'
    // 2.25.0_12312025 - must edit manually below (line 28)!!!
    container 'quay.io/jvhagey/mlst@sha256:a67904d356118f9c163d26000d4d78cc449e3205145f87be28726869d67602f7'

    input:
    tuple val(meta), path(fasta), path(taxonomy)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path("versions.yml")          , emit: versions

    script:
        //set up for terra
    if (params.terra==false) {
        terra_path = ""
        terra_exit = ""
    } else if (params.terra==true) {
        terra_path = "PATH=/opt/conda/envs/mlst/bin:\$PATH"
        terra_exit = """PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/mlst/bin:||')" """
    } else {
        error "Please set params.terra to either \"true\" or \"false\""
    }
    // helps set correct paths to get database version being used
    def terra = params.terra ? "true" : "false"
    //define variables
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // mlst is suppose to allow gz and non-gz, but when run in the container (outside of the pipeline) it doesn't work. Also, doesn't work on terra so adding unzip step
    def container = task.container.toString() - "quay.io/jvhagey/mlst@"
    def mlst_version = "2.25.0_20251231"
    def mlst_version_clean = mlst_version.split("_")[0]
    """
    #adding mlst path for running mlst on terra
    $terra_path

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

    # Address all the cases where we need to run mlst twice to get both schemes
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
    elif [[ \$scheme == "aparagallinarum_Ghanem" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme aparagallinarum_Guo --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "aparagallinarum_Guo" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme aparagallinarum_Ghanem --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "efaecium" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme efaecium_Bezdicek --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "efaecium_Bezdicek" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme efaecium --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "leptospira" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme leptospira_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        mlst --scheme leptospira_3 --threads $task.cpus \$unzipped_fasta > ${prefix}_3.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv ${prefix}_3.tsv> ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "leptospira_2" ]]; then
        mv ${prefix}.tsv ${prefix}_2.tsv
        mlst --scheme leptospira --threads $task.cpus \$unzipped_fasta > ${prefix}_1.tsv
        mlst --scheme leptospira_3 --threads $task.cpus \$unzipped_fasta > ${prefix}_3.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv ${prefix}_3.tsv> ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "leptospira_3" ]]; then
        mv ${prefix}.tsv ${prefix}_3.tsv
        mlst --scheme leptospira_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        mlst --scheme leptospira --threads $task.cpus \$unzipped_fasta > ${prefix}_1.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv ${prefix}_3.tsv> ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "salmonella_Oxford" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme salmonella_Achtman --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "salmonella_Achtman" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme salmonella_Oxford --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "mgallisepticum_Ghanem" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme mgallisepticum_Beko --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "mgallisepticum_Beko" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme mgallisepticum_Ghanem --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "pmultocida_multihost" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme pmultocida_rirdc --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "pmultocida_rirdc" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme pmultocida_multihost --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "mbovis" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme mbovis_legacy --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "mbovis_legacy" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme mbovis --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "smutans_Do" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme smutans_Kakano --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "smutans_Kakano" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme smutans_Do --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "sthermophilus" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme sthermophilus_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "sthermophilus_2" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme sthermophilus --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "tpallidum_Grillova" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme tpallidum_Pla-Diaz --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \$scheme == "tpallidum_Pla-Diaz" ]]; then
        mv ${prefix}.tsv ${prefix}_1.tsv
        mlst --scheme tpallidum_Grillova --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
        cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
        rm ${prefix}_*.tsv
    elif [[ \${genus,,} == "mycobacterium" ]] || [[ \${genus,,} == "mycobacteroides" ]] ; then
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
        elif [[ \$scheme == "citrobacter" ]]; then
            mv ${prefix}.tsv ${prefix}.OLD-tsv
            mlst --scheme ecoli --threads $task.cpus \$unzipped_fasta > ${prefix}_1.tsv
            mlst --scheme ecoli_2 --threads $task.cpus \$unzipped_fasta > ${prefix}_2.tsv
            cat ${prefix}_1.tsv ${prefix}_2.tsv > ${prefix}.tsv
            rm ${prefix}_*.tsv
        elif [[ \$scheme == "salmonella_Achtman" ]] || [[ \$scheme == "salmonella_Oxford" ]]; then
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
    sed -i '1i source_file  Database  ST  locus_1 locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 locus_8 locus_9  locus_10' ${prefix}.tsv

    #handling to get database version being used
    if [[ $terra == false ]]; then
        db_version=\$(head -n 1 "/mlst-${mlst_version_clean}/db/db_version" | awk '{\$1=\$1; print}')
    else
        db_version=\$(cat /opt/conda/envs/phoenix/db/db_version | date -f - +%Y-%m-%d )
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( echo \$(mlst --version 2>&1) | sed 's/mlst //' )
        mlst_db: \$db_version
        mlst_container: ${container}
    END_VERSIONS

    #revert path back to main envs for running on terra
    $terra_exit
    """
}