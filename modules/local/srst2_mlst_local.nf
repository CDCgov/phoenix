process SRST2_MLST {
	tag "${meta.id}"
	label 'process_medium'
	// 0.2.0_patched
	container 'quay.io/jvhagey/srst2@sha256:d4a68baf84c8818b59f334989ccbeea044baf14a79aaf1bd95a1f24f69d0dc5b'

	input:
	tuple val(meta), path(fastqs), path(getmlstout), path(alleles), path(profiles), val(status)

	output:
	//tuple val(meta), path("*_mlst_*_results.txt")				, optional:true, emit: mlst_results
	tuple val(meta), path("*_srst2.mlst")						, optional:true, emit: mlst_results
	tuple val(meta), path("*_srst2_temp.mlst") 					, optional:true, emit: mlst_results_temp
	tuple val(meta), path("*.pileup")							, optional:true, emit: pileup
	tuple val(meta), path("*.sorted.bam")						, optional:true, emit: sorted_bam
	tuple val(meta), path("*_srst2_status.txt")					, optional:true, emit: empty_checker
	path "versions.yml"											,				emit: versions

	when:
	(task.ext.when == null || task.ext.when)

	script:
	// set up terra variables
	if (params.terra==false) {
		terra = ""
		terra_exit = ""
	} else if (params.terra==true) {
		terra = """export PYTHONPATH=/opt/conda/envs/srst2/lib/python2.7/site-packages/
		PATH=/opt/conda/envs/srst2/bin:\$PATH
		"""
		terra_exit = """export PYTHONPATH=/opt/conda/envs/phoenix/lib/python3.7/site-packages/
		PATH="\$(printf '%s\\n' "\$PATH" | sed 's|/opt/conda/envs/srst2/bin:||')"
		"""
	} else {
		error "Please set params.terra to either \"true\" or \"false\""
	}
	// define variables
	def args = task.ext.args ?: ""
	def prefix = task.ext.prefix ?: "${meta.id}"
	def read_s = meta.single_end ? "--input_se ${fastqs}" : "--input_pe ${fastqs[0]} ${fastqs[1]}"
	def container = task.container.toString() - "quay.io/jvhagey/srst2@"
	
	"""
	#adding python path for running srst2 on terra
	$terra

	echo "STATUS-IN: ${status[0]}"
	srst2_ran="False"
	if [[ "${status[0]}" = "False" ]]; then
		scheme_count=1
		for getout in $getmlstout
		do
			no_match="False"
			echo "\${getout}"
			line="\$(tail -n1 \${getout})"
			if [[ "\${line}" = "DB:No match found"* ]] || [[ "\${line}" = "DB:Server down"* ]]; then
				no_match="True"
				mlst_db="No_match_found"
			else
				mlst_db=\$(echo "\${line}" | cut -f1 | cut -d':' -f2)
				# Put check much closer to where it is needed. Headers and labels dont match due to db inclusion. Removed to make srst2 happy

				mlst_delimiter=\$(echo "\${line}" | cut -f3 | cut -d':' -f2 | cut -d"'" -f2)

				echo "Test: \${mlst_db} \${mlst_db}_profiles.csv \${mlst_delimiter}"

				srst2 ${read_s} \\
					--threads $task.cpus \\
					--output \${scheme_count}_${prefix} \\
					--mlst_db \${mlst_db}_temp.fasta \\
					--mlst_definitions \${mlst_db}_profiles_temp.csv \\
					--mlst_delimiter \${mlst_delimiter} \\
					$args
			fi

			if [[ -f \${scheme_count}_${prefix}__mlst__\${mlst_db}_temp__results.txt ]]; then
				lines_in_result_file=\$(cat \${scheme_count}_${prefix}__mlst__\${mlst_db}_temp__results.txt | wc -l)
			else
				lines_in_result_file=0
				echo "No srst2 result file exists"
			fi

			header="Sample	database	ST	mismatches	uncertainty	depth	maxMAF	locus_1	locus_2	locus_3	locus_4	locus_5	locus_6	locus_7	locus_8	locus_9	locus_10"
			if [[ "\${scheme_count}" -eq 1 ]]; then
				echo "\${header}" > ${prefix}_srst2.mlst
			fi

			header_list=""
			trailer_list=""
			if [[ "\${no_match}" = "True" ]]; then
				tax_with_no_scheme=\$(echo "\${line}" | cut -d'(' -f2 | cut -d')' -f1)
				echo "${prefix}	-	-	-	-	-	-" >> "${prefix}_srst2.mlst"
			elif [[ "\${lines_in_result_file}" -eq 1 ]]; then
				echo "Not enough was found to even make a guess"
#				echo "${prefix}	No match found for \${mlst_db}	-	-	-	-	-" >> "${prefix}_srst2.mlst"
				echo "${prefix}	-	-	-	-	-	-" >> "${prefix}_srst2.mlst"
			else
				raw_header="\$(head -n1 \${scheme_count}_${prefix}*.txt)"
				# Account for the cases where multi-databases require extra genes ID's during processing, but remove here.
				# Current known list is only populated by Abaumannii
				to_remove_prefixes=('Pas_' 'Ox_')
				trimmed_header="\${raw_header}"
				for prefix in \${to_remove_prefixes[@]}; do
					trimmed_header="\${trimmed_header//\${prefix}/}"
					echo "\${prefix}, \${trimmed_header}"
				done
				raw_trailer="\$(tail -n1 \${scheme_count}_${prefix}*.txt)"
				formatted_trailer="${prefix}	\${mlst_db}"
				IFS=\$'\t' read -r -a trailer_list <<< "\$raw_trailer"
				IFS=\$'\t' read -r -a header_list <<< "\$trimmed_header"
				header_length="\${#header_list[@]}"
				ST_index=1
				mismatch_index=\$(( header_length - 4 ))
				uncertainty_index=\$(( header_length - 3 ))
				depth_index=\$(( header_length - 2 ))
				maxMAF_index=\$(( header_length - 1))
				genes_start_index=2
				genes_end_index=\$(( header_length - 5 ))
				# Original Index will be something as follows for a 7 gene scheme
				#
				# Sample  ST	  adk	 atpA	dxr	 glyA	recA	sodA	tpi	 mismatches	  uncertainty	 depth   maxMAF
				# 0	   1	   2	   3	   4	   5	   6	   7	   8	   9				10			 11	  12
				#
				# Expected Index will be as follows for the same isolate to accomomdate splunk ingestion of varying gene counts
				#
				# Sample  database  ST  mismatches  uncertainty depth maxMAF  locus_1 locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 locus_8 locus_9 locus_10
				# 0	 1   2   3		   4		   5	 6	   7	   8	   9	   10	  11	  12		13	14	  15	  16

				echo "\${#header_list[@]} --- \${header_list[@]} --- \${#trailer_list[@]} --- \${trailer_list[@]}"
				formatted_trailer="\${formatted_trailer}	\${trailer_list[\${ST_index}]}"
				formatted_trailer="\${formatted_trailer}	\${trailer_list[\${mismatch_index}]}"
				formatted_trailer="\${formatted_trailer}	\${trailer_list[\${uncertainty_index}]}"
				formatted_trailer="\${formatted_trailer}	\${trailer_list[\${depth_index}]}"
				formatted_trailer="\${formatted_trailer}	\${trailer_list[\${maxMAF_index}]}"

				#for index in {\$genes_start_index..\$genes_end_index}
				for (( index=\${genes_start_index} ; index <= \${genes_end_index} ; index++ ));
				do
					echo "\${index} -- \${header_list[\${index}]} -- \${trailer_list[\${index}]}"
					formatted_trailer="\${formatted_trailer}	\${header_list[\${index}]}(\${trailer_list[\${index}]})"
				done
				echo "\${formatted_trailer}" >> ${prefix}_srst2.mlst
			fi

			scheme_count=\$(( scheme_count + 1 ))
			#mv profiles_csv \${mlst_db}_profiles_csv
			cp ${prefix}_srst2.mlst ${prefix}_srst2_temp.mlst
			srst2_ran="True"
		done
	else
		echo "DONT USE" > ${prefix}_srst2_temp.mlst
	fi

	echo "\${srst2_ran}" > "${prefix}_srst2_status.txt"

	cat <<-END_VERSIONS > versions.yml
	"${task.process}":
	    srst2: \$(echo \$(srst2 --version 2>&1) | sed 's/srst2 //' )
	    srst2_commit_patched: 73f885f55c748644412ccbaacecf12a771d0cae9
	    srst2_container: ${container}
	END_VERSIONS

	#revert python path back to main envs for running on terra
	$terra_exit
	"""
}
	
