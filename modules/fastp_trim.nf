/*

Fastp trim module

*/


process FASTP_TRIM {

  tag {row.sm}

  publishDir "${params.output}/metrics", mode: "copy", pattern: "*_fastp.html"

  label "large"

  input:
    tuple val(row), path(fq1), path(fq2)

  output:
    tuple val(row), path(fq1), path(fq2), emit: fastq_trimmed
    path "*_fastp.json", emit: fastp_json
    path "*_fastp.html"

  """
    fastp -i ${fq1} -I ${fq2} \\
      -o ${row.sm}_1P.fq.gz \\
      -O ${row.sm}_2P.fq.gz \\
      --length_required 20 \\
      -j ${row.sm}_fastp.json \\
      -h ${row.sm}_fastp.html
  """
}
