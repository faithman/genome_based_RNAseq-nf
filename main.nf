#!/usr/bin/env nextflow

/*
    gbRNA-seq Pipeline
    Chengdu Research Base of Giant Panda Breeding
    Authors:
    - Ye Wang <yewangfaith@gmail.com>
*/
nextflow.enable.dsl=2

// This pipeline requires NXF_VER 20.01.0-rc1 or later

/*
    Params
*/


date = new Date().format( 'yyyyMMdd' )
params.debug =  false
params.email = "yewang_faith@hotmail.com"

// Check that reference exists
params.genome_dir = ""
params.gtf = ""
params.cpus = ""
params.memory = ""

//reference = file(params.reference, checkIfExists: true)

// Debug
if (params.debug.toString() == "true") {
    params.output = "debug-${date}"
    params.sample_sheet = "${workflow.projectDir}/test_data/test_sample_sheet.tsv"
} else {
    // The strain sheet that used for 'production' is located in the root of the git repo
    params.output = "standard-${date}"
    params.sample_sheet = "${workflow.projectDir}/sample_sheet.tsv"
}


def log_summary() {

    out =  '''

    ▐   ▗▄▄ ▗▖ ▖ ▗▖     ▗▖ ▖▗▄▄▖
 ▄▄ ▐▄▖ ▐ ▝▌▐▚ ▌ ▐▌     ▐▚ ▌▐
▐▘▜ ▐▘▜ ▐▄▄▘▐▐▖▌ ▌▐     ▐▐▖▌▐▄▄▖
▐ ▐ ▐ ▐ ▐ ▝▖▐ ▌▌ ▙▟  ▀▘ ▐ ▌▌▐
▝▙▜ ▐▙▛ ▐  ▘▐ ▐▌▐  ▌    ▐ ▐▌▐
 ▖▐
 ▝▘
                                              
'''

out += """
    parameters               description                 Set/Default
    ==========               ===========                 ========================
    Debug                    Debug or not                ${params.debug}
    output                   Release Directory           ${params.output}
    sample_sheet             sample sheet                ${params.sample_sheet}
    genome_dir               Genome and index            ${params.genome_dir}
    username                                             ${"whoami".execute().in.text}

    Nextflow Run
    ---------------
    ${workflow.commandLine}
    run name                                             ${workflow.runName}
    scriptID                                             ${workflow.scriptId}
    git commit                                           ${workflow.commitId}
    container                                            ${workflow.container}
---
"""
out
}

log.info(log_summary())


// Includes processes from modules
include { FASTP_TRIM } from './modules/fastp_trim.nf' params(params)
include { STAR } from './modules/star.nf' params(params)
include { PICARD_MARKDUPLICATES } from './modules/picard_markdup.nf' params(params)
include { SAMBAMBA_INDEX } from './modules/sambamba_index.nf' params(params)
include { SAMTOOL_STATS } from './modules/samtool_stats.nf' params(params)
include { STRINGTIE_QUANT } from './modules/stringtie.nf' params(params)


if (workflow.profile == "") {
    println "Must set -profile: local, pbs, slurm"
    exit 1
}

// Read sample sheet into channel
sample_sheet = Channel 
    .fromPath(params.sample_sheet, checkIfExists: true)
    .ifEmpty {exit 1, "sample sheet not found"}
    .splitCsv(header:true, sep: "\t")
    .map { row -> row.fq1 = row.fq1;row}
    .map { row -> row.fq2 = row.fq2;row}
    .map { row -> [row, file(row.fq1), file(row.fq2)]}
    //.view()

    
// Workflow
workflow {

    // Generate a summary of the current run
    summary(Channel.from("run"))

    // Trim raw fastq files
    FASTP_TRIM( sample_sheet )

    // Perform mapping
    STAR( FASTP_TRIM.out.fastq_trimmed, params.genome_dir )

    // Markduplicates, for pair-end reads, it definitely need to mark duplicates 
    PICARD_MARKDUPLICATES( STAR.out.bam_sorted )

    // Indexing bam
    SAMBAMBA_INDEX( PICARD_MARKDUPLICATES.out.bam_marked )

    // Bam stats
    SAMTOOL_STATS( PICARD_MARKDUPLICATES.out.bam_marked.join(SAMBAMBA_INDEX.out.bai, by: [0]) )

    // Quant
    STRINGTIE_QUANT( PICARD_MARKDUPLICATES.out.bam_marked.join(SAMBAMBA_INDEX.out.bai, by: [0]), params.gtf )
}


process summary {

    // This process is for a summary of current run.
    
    executor 'local'

    publishDir "${params.output}", mode: 'copy'
    
    input:
        val(run)

    output:
        path("sample_sheet.tsv")
        path("summary.txt")

    """
        echo '''${log_summary()}''' > summary.txt
        cat ${params.sample_sheet} > sample_sheet.tsv
    """

}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: params.email, subject: 'My pipeline execution', body: msg)
}
