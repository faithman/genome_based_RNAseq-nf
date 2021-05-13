/*

Picard markduplicates module

*/

process PICARD_MARKDUPLICATES {
	
	tag {sm}

	publishDir "${params.output}/bams", mode: "copy", pattern: "*.bam"
	publishDir "${params.output}/metrics", mode: "copy", pattern: "*.metrics.txt"

	label "large"

	input:
	tuple val(sm), path(bam)

	output:
	tuple val(sm), path("*.bam"), emit: bam_marked
	tuple val(sm), path("*.metrics.txt"), emit: metrics

	"""
	picard \\
		-Xmx8g \\
		MarkDuplicates \\
		INPUT=${bam} \\
		OUTPUT=${sm}.bam \\
		METRICS_FILE=${sm}.MarkDuplicates.metrics.txt
	"""
}
