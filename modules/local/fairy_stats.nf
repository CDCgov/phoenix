process FAIRY_STATS {
	tag "$meta.id"
	label 'process_medium'
	container 'quay.io/jvhagey/phoenix:base_v1.1.0'

	input:
	tuple val(meta), path(reads)

	output:
	tuple val(meta), path('*_stats.txt'),optional:true,                     emit: raw_stats
	tuple val(meta), path('*_raw_read_counts.txt'),optional:true,           emit: combined_raw_stats
	path("versions.yml"),optional:true,                                     emit: versions
	tuple val(meta), path('*_result.txt'),optional:true,                    emit: outcome
	tuple val(meta), path('*.fastq.gz'),optional:true,                      emit: reads
	tuple val(meta),path('*_prdresult.txt'),optional:true,                  emit: compr_rds
	path('*_summaryline_failure.tsv'),optional:true,                        emit: summary_line
	path('*.synopsis'),optional:true,                        emit: synopsis

	script: // This script is bundled with the pipeline, in cdcgov/phoenix/bin/
	def prefix = task.ext.prefix ?: "${meta.id}"
	def num1 = "${reads[0]}".minus(".fastq.gz")
	def num2 = "${reads[1]}".minus(".fastq.gz")
	"""
	set +e
fairy_proc.sh ${reads[0]} 
fairy_proc.sh ${reads[1]}
if grep "PASS" ${prefix}_results.txt; then
	q30.py ${reads[0]} > ${prefix}_R1_stats.txt
	q30.py ${reads[1]} > ${prefix}_R2_stats.txt
fi

if [ -f ${prefix}_R1_stats.txt -a -f ${prefix}_R2_stats.txt ]; then
		create_raw_stats_output.py -n ${prefix} -r1 ${prefix}_R1_stats.txt -r2 ${prefix}_R2_stats.txt
		fairy.py -r ${prefix}_raw_read_counts.txt
fi

if [ ! -f ${prefix}_summaryline_failure.tsv ] ; then
	mv ${reads[0]} ${num1}_C.fastq.gz
	mv ${reads[1]} ${num2}_C.fastq.gz
else
	if [ -f ${prefix}_raw_read_counts.txt ]; then
		#create synopsis file
		fairy_countdiff.sh ${prefix}_raw_read_counts.txt
	fi

	#delete files that would enter the channel
	rm ${reads[0]}
	rm ${reads[1]}
	mv ${prefix}_result.txt ${prefix}_prdresult.txt
	rm ${prefix}_prdresult.txt
fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    python: \$(python --version | sed 's/Python //g')
END_VERSIONS
	"""
}
