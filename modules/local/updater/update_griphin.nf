process UPDATE_GRIPHIN {
    tag "${full_project_id.toString().split('/')[-1].replace("]","")}"
    label 'process_low'
    stageInMode 'copy' // you need this or openpyxl complains that excel files aren't excel files. 
    container 'quay.io/jvhagey/phoenix@sha256:ba44273acc600b36348b96e76f71fbbdb9557bb12ce9b8b37787c3ef2b7d622f'

    input:
    path(griphins_excel)
    //path(griphins_tsv)
    val(full_project_id) // output directory used as prefix for the summary file
    path(valid_samplesheet_file)
    val(coverage)
    path(bldb)
    val(remove_dups_var)
    // This may need edited as all situations may not have been accounted for but if a list is used then a 'final project' name needs to be used. Typically, the outdir or griphin_out folder location, for species specific pipelines.
    // If only 2 are given then the name will just stay as the original name of the file.
    val(outdir_location) // Is just the name of where the file will go, different from the meta.project_id

    output:
    path("*_GRiPHin_Summary.xlsx"), emit: griphin_report
    path("*_GRiPHin_Summary.tsv"),  emit: griphin_tsv_report
    path("versions.yml"),           emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    def ica = params.ica ? "python ${params.bin_dir}" : ""
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def project_path = file(full_project_id).parent // only removed the last dir, but gives the full path otherwise
    def project_id = full_project_id.toString().split('/')[-1] // just the last dir name 
    // passing in either two files or a list of files
    def griphin_input = griphins_excel.size() == 2 ? "-g1 ${griphins_excel[0]} -g2 ${griphins_excel[1]}" : "--griphin_list"
    // file name handling
    def out_name = (params.mode_upper == "UPDATE_PHOENIX" || griphins_excel.size() > 2) ? "${outdir_location}_GRiPHin_Summary" : "${project_id}_GRiPHin_Summary"
    /*if (griphins_excel.size() == 2) { // if there are only two files then we aren't dealing with multi_dir situation 
        // Case where only two files are passed
        griphin_input = "-g1 ${griphins_excel[0]} -g2 ${griphins_excel[1]}"
        out_name = "${project_id}_GRiPHin_Summary" //--> not working for updater, --input with 2 dirs, gives full path name, which I ends up published wrong
    } else if (griphins_excel.size() > 2) {
        // Case where griphins_excel contains many
        griphin_input = "--griphin_list"
        out_name = "${outdir_location}_GRiPHin_Summary"
    } else {
        error "The input griphins_excel must contain at least two files to combine."
    }*/
    def valid_samplesheet = valid_samplesheet_file ? "--samplesheet ${valid_samplesheet_file}" : ""
    def remove_dups = remove_dups_var ? "--remove_dups" : "" // When we are running species specific and updater pipelines we need to remove duplications from the original griphin reports so we can update with new data
    """
    ${ica}combine_GRiPHins.py ${griphin_input} \
        --output ${out_name} \
        --coverage ${coverage} \
        --parent_folder ${project_path} \
        --bldb ${bldb} \
        ${remove_dups} \
        ${valid_samplesheet}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        combine_GRiPHins.py: \$(${ica}combine_GRiPHins.py --version)
        phoenix_base_container_tag: ${container_version}
        phoenix_base_container: ${container}
    END_VERSIONS
    """
}