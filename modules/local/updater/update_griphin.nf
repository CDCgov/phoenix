process UPDATE_GRIPHIN {
    tag "${full_project_id.toString().split('/')[-1].replace("]","")}"
    label 'process_low'
    stageInMode 'copy' // you need this or openpyxl complains that excel files aren't excel files. 
    container 'quay.io/jvhagey/phoenix@sha256:2122c46783447f2f04f83bf3aaa076a99129cdd69d4ee462bdbc804ef66aa367'

    input:
    path(griphins_excel)
    //path(griphins_tsv)
    path(full_project_id) // output directory used as prefix for the summary file
    path(valid_samplesheet_file)
    val(coverage)
    path(bldb)
    val(remove_dups_var)
    // This may need edited as all situations may not have been accounted for but if a list is used then a 'final project' name needs to be used. Typically the outdir or griphin_out folder location.
    // If only 2 are given then the name will just stay as the original name of the file.
    val(outdir_location) // Is just the name of where the file will go, different from the meta.project_id

    output:
    path("*_GRiPHin_Summary.xlsx"), emit: griphin_report
    path("*_GRiPHin_Summary.tsv"),  emit: griphin_tsv_report
    path("versions.yml"),           emit: versions

    script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
    // Adding if/else for if running on ICA it is a requirement to state where the script is, however, this causes CLI users to not run the pipeline from any directory.
    if (params.ica==false) { ica = "" } 
    else if (params.ica==true) { ica = "python ${params.bin_dir}" }
    else { error "Please set params.ica to either \"true\" if running on ICA or \"false\" for all other methods." }
    // define variables
    def container_version = "base_v2.2.0"
    def container = task.container.toString() - "quay.io/jvhagey/phoenix@"
    def project_id = full_project_id.toString().split('/')[-1].replace("]","")
    if (griphins_excel.size() == 2) {
        // Case where only two files are passed
        griphin_input = "-g1 ${griphins_excel[0]} -g2 ${griphins_excel[1]}"
        out_name = "${project_id}_GRiPHin_Summary"
    } else {
        // Case where griphins_excel contains many
        griphin_input = "--griphin_list"
        out_name = "${outdir_location}_GRiPHin_Summary"
    }
    def valid_samplesheet = valid_samplesheet_file ? "--samplesheet ${valid_samplesheet_file}" : ""
    def remove_dups = remove_dups_var ? "--remove_dups" : "" // When we are running species specific and updater pipelines we need to remove duplications from the original griphin reports so we can update with new data
    """
    ${ica}combine_GRiPHins.py ${griphin_input} \
        --output ${out_name} \
        --coverage ${coverage} \
        --parent_folder ${project_id} \
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