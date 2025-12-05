process GRIPHIN {
    label 'process_low'
    // base_v2.3.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:b8e3d7852e5f5b918e9469c87bfd8a539e4caa18ebb134fd3122273f1f412b05'

    input:
    path(db)
    path(original_samplesheet)
    val(metas)
    path(griphin_files) // path to the GRiPHin files, which includes the staged directory
    path(outdir) //output directory used as prefix for the summary file
    val(phx_version)
    val(coverage)
    val(phx_mode)
    val(shigapass_detected)
    val(centar_detected)
    path(bldb)
    val(filter_var) // needed for species specific mode
    val(dont_publish)
    path(blind_list) // -b (mostly needed for PhyloPHoenix runs)

    output:
    tuple path("full_path_file.txt"), path("*_GRiPHin*.xlsx"),                         emit: griphin_report
    tuple path("full_path_file.txt"), path("*_GRiPHin*.xlsx"), path("*_GRiPHin*.tsv"), emit: griphins
    path("*_GRiPHin*.tsv"),                                                            emit: griphin_tsv_report
    path("Directory_samplesheet.csv"), optional:true,                                  emit: converted_samplesheet //the only time this isn't made is with --centar with --samplesheet
    path("versions.yml"),                                                              emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def blind_names   = blind_list ? "--blind_list ${blind_list}" : ""
    def phoenix_mode = phx_mode ? "" : "--phoenix"
    def scaffolds = (params.mode_upper == "SCAFFOLDS" || params.mode_upper == "CDC_SCAFFOLDS") ? "--scaffolds" : "" 
    def shigapass = shigapass_detected ? "--shigapass" : ""
    def centar = centar_detected ? "--centar" : ""
    def updater = (params.mode_upper == "UPDATE_PHOENIX") ? "--updater" : "" 
    //def samplesheet_command = (centar_detected && original_samplesheet) ? "--samplesheet ${original_samplesheet}" : ""
    def filter = filter_var ? "--filter_samples" : ""
    def output_prefix = ((dont_publish == true) || (params.mode_upper == "CENTAR" && params.indir == null)) ? "${outdir}_GRiPHin" : "${outdir}_GRiPHin_Summary" 
    def container_version = "base_v2.3.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    def prefix = task.ext.prefix ?: "GRiPHin"
    def stage_files = [
        metas.collect { "mkdir -p ${prefix}/${it.id}" },
        metas.collect { "mv ${it.filenames.join(' ')} ${prefix}/${it.id}" }
    ].flatten().join(" && ")
    """
    ${stage_files}

    # Get full path to outdir for parent folder and data location columns in griphin
    full_path=\$(readlink -f ${outdir})

    # Save the full_path to a file (this file will be captured in the output block)
    echo \$full_path > full_path_file.txt

    ${ica}GRiPHin.py -d \$full_path -a $db --output ${output_prefix} --bldb ${bldb} --phx_version ${phx_version} ${filter}\
        --coverage ${coverage} ${phoenix_mode} ${shigapass} ${centar} ${scaffolds} --samplesheet ${original_samplesheet} ${updater} ${blind_names}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        griphin.py: \$(${ica}GRiPHin.py --version)
        bldb_creation_date: \$(echo ${bldb} | sed 's/[^0-9]//g')
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}