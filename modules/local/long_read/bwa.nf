process BWA {
  tag "${meta.id}"
  label 'process_medium'
  container 'staphb/bwa:0.7.17'
  //sha256:5352be51d07f011974dac4c0db8800731360d435dfd0a6258d1cd1c877166bd1
  errorStrategy 'ignore'

  input:
  tuple val(meta), file(fasta), file(reads)

  output:
  tuple val(meta), file("${meta.id}_{1,2}.sam"), emit: sam
  path "versions.yml",                           emit: versions

  script:
  """
    bwa index $fasta
    bwa mem -t 16 -a $fasta ${meta.id}_1.trim.fastq.gz > ${meta.id}_1.sam
    bwa mem -t 16 -a $fasta ${meta.id}_2.trim.fastq.gz > ${meta.id}_2.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
  """
}
