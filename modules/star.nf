/*

STAR mapping module

*/


process STAR {

  tag {row.sm}

  //publishDir "${params.output}/bams", mode: "copy"

  //label "large"

  cpus 4 
  memory 8

  input:
    tuple val(row), path(fq1), path(fq2)
    path genome_dir

  output:
    tuple val(row.sm), path("*Aligned.out.bam"), optional:true, emit: bam
    tuple val(row.sm), path("*sortedByCoord.out.bam"), optional:true, emit: bam_sorted
    tuple val(row.sm), path("*fastq.gz"), optional:true, emit: fastq


  """
    STAR \\
      --runThreadN ${task.cpus} \\
      --genomeDir ${genome_dir} \\
      --readFilesIn ${fq1} ${fq2} \\
      --outSAMtype BAM SortedByCoordinate \\
      --outReadsUnmapped Fastx \\
      --readFilesCommand zcat \\
      --outSAMattrRGline ID:${row.group} SM:${row.sm} \\
      --outFileNamePrefix ${row.sm}_ 


    if [ -f ${row.sm}.Unmapped.out.mate1 ]; then
        mv ${row.sm}.Unmapped.out.mate1 ${row.sm}.unmapped_1.fastq
        bgzip ${row.sm}.unmapped_1.fastq
    fi
    if [ -f ${row.sm}.Unmapped.out.mate2 ]; then
        mv ${row.sm}.Unmapped.out.mate2 ${row.sm}.unmapped_2.fastq
        bgzip ${row.sm}.unmapped_2.fastq
    fi
  """
}
