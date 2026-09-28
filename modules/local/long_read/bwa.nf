process BWA {
  tag "${meta.id}"
  label 'process_medium'
  container 'staphb/bwa@sha256:e4f2bd6ba48ad1923f2edec641a59f2ba0f26b59805fb9812d780996ba5fa8df'

  input:
  tuple val(meta), file(fasta), file(reads)

  output:
  tuple val(meta), file("${meta.id}_{1,2}.sam"), emit: sam
  path "versions.yml",                           emit: versions

  script:
  def container = task.container.toString() - "staphb/bwa@"
  """
    bwa index $fasta
    bwa mem -t 16 -a $fasta ${meta.id}_1.trim.fastq.gz > ${meta.id}_1.sam
    bwa mem -t 16 -a $fasta ${meta.id}_2.trim.fastq.gz > ${meta.id}_2.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        bwa_container: ${container}
    END_VERSIONS
  """
}
