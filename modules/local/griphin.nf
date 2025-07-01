process GRIPHIN {
    label 'process_low'
    // base_v2.2.0 - MUST manually change below (line 36)!!!
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    path(summary_line_files)
    path(original_samplesheet)
    path(db)
    tuple path(outdir), path(dir2) // output directory used as prefix for the summary file, dir2 is the location of original input folder (need for update_phoenix when multi dirs are given in --input )
    val(phx_version)
    val(coverage)
    val(entry)
    val(shigapass_detected)
    val(centar_detected)
    path(bldb)
    val(filter_var) // needed for species specific entry points

    output:
    tuple path("full_path_file.txt"), path("*_GRiPHin*.xlsx"),                         emit: griphin_report
    tuple path("full_path_file.txt"), path("*_GRiPHin*.xlsx"), path("*_GRiPHin*.tsv"), emit: griphins
    path("*_GRiPHin*.tsv"),                                                            emit: griphin_tsv_report
    path("Directory_samplesheet.csv"),          optional:true,                         emit: converted_samplesheet //the only time this isn't made is with --centar with --samplesheet
    path("versions.yml"),                                                              emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def phoenix = entry ? "--phoenix" : ""
    def scaffolds = (params.pipeline_upper == "SCAFFOLDS" || params.pipeline_upper == "CDC_SCAFFOLDS") ? "--scaffolds" : "" 
    def shigapass = shigapass_detected ? "--shigapass" : ""
    def centar = centar_detected ? "--centar" : ""
    def updater = (params.pipeline_upper == "UPDATE_PHOENIX") ? "--updater" : "" 
    //def samplesheet_command = (centar_detected && original_samplesheet) ? "--samplesheet ${original_samplesheet}" : ""
    def filter = filter_var ? "--filter_samples" : ""
    def output_prefix = ((params.pipeline_upper == "UPDATE_PHOENIX" && params.indir == null) || params.pipeline_upper == "CENTAR") ? "${outdir}_GRiPHin" : "${outdir}_GRiPHin_Summary" 
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix:"
    """
    full_path=\$(readlink -f ${outdir})

    # Save the full_path to a file (this file will be captured in the output block)
    echo \$full_path > full_path_file.txt

    ${ica}GRiPHin.py -d \$full_path -a $db --output ${output_prefix} --bldb ${bldb} --phx_version ${phx_version} ${filter}\
        --coverage ${coverage} ${phoenix} ${shigapass} ${centar} ${scaffolds} --samplesheet ${original_samplesheet} ${updater}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        griphin.py: \$(${ica}GRiPHin.py --version)
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}