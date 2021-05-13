/*

Samtool stats module

*/

process SAMTOOL_STATS {

    tag {sm}

    publishDir "${params.output}/stats", mode: "copy", pattern: "*.stats"

    label "large"

    input:
    tuple val(sm), path(bam), path(bai)
    
    output:
    tuple val(sm), path("*.stats"), emit: stats


    """
    samtools stats ${bam} > ${bam}.stats
    """
}