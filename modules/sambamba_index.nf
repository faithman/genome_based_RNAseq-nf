/*

Sambamba index module

*/

process SAMBAMBA_INDEX {
	
	tag {sm}

	publishDir "${params.output}/bams", mode: "copy", pattern: "*.bai"

	label "large"

	input:
	tuple val(sm), path(bam)

	output:
	tuple val(sm), path("*.bai"), emit: bai

	"""
	sambamba \\
		index ${bam}
	"""
}