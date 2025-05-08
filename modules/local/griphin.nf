process GRIPHIN {
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    path(summary_line_files)
    path(original_samplesheet)
    path(db)
    path(outdir) // output directory used as prefix for the summary file
    val(phx_version)
    val(coverage)
    val(entry)
    val(scaffolds_entry)
    val(update) //should be true to for updater and species specific entry points
    val(shigapass_detected)
    val(centar_detected)
    path(bldb)
    val(filter_var) // needed for species specific entry points

    output:
    tuple path("full_path_file.txt"), path("*_GRiPHin*.xlsx"), emit: griphin_report
    path("*_GRiPHin*.tsv"),                                    emit: griphin_tsv_report
    path("Directory_samplesheet.csv"),          optional:true, emit: converted_samplesheet //the only time this isn't made is with --centar with --samplesheet
    path("versions.yml"),                                      emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def phoenix = entry ? "--phoenix" : ""
    def scaffolds = scaffolds_entry ? "--scaffolds" : ""
    def shigapass = shigapass_detected ? "--shigapass" : ""
    def centar = centar_detected ? "--centar" : ""
    //def samplesheet_command = (centar_detected && original_samplesheet) ? "--samplesheet ${original_samplesheet}" : ""
    def filter = filter_var ? "--filter_samples" : ""
    def output_prefix = update ? "${outdir}_GRiPHin" : "${outdir}_GRiPHin_Summary"
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    full_path=\$(readlink -f ${outdir})

    # Save the full_path to a file (this file will be captured in the output block)
    echo \$full_path > full_path_file.txt

    ${ica}GRiPHin.py -d \$full_path -a $db --output ${output_prefix} --bldb ${bldb} --phx_version ${phx_version} ${filter}\
        --coverage ${coverage} ${phoenix} ${shigapass} ${centar} ${scaffolds} --samplesheet ${original_samplesheet}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        griphin.py: \$(${ica}GRiPHin.py --version)
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}