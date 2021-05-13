/*

Stringtie quant module

*/


process STRINGTIE_QUANT {

  tag {sm}

  publishDir "${params.output}/quant", mode: "copy"

  label "large"

  input:
    tuple val(sm), path(bam), path(bai)
    path gtf

  output:
    tuple val(sm), path("*.coverage.gtf")   , emit: coverage_gtf
    tuple val(sm), path("*.transcripts.gtf"), emit: transcript_gtf
    tuple val(sm), path("*.txt")            , emit: abundance
    tuple val(sm), path("*.ballgown")       , emit: ballgown


  """
    stringtie \\
        $bam \\
        -G ${gtf} \\
        -o ${sm}.transcripts.gtf \\
        -A ${sm}.gene_abundance.txt \\
        -C ${sm}.coverage.gtf \\
        -b ${sm}.ballgown 
  """
}
