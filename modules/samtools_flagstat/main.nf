#!/usr/bin/env nextflow

process SAMTOOLS_FLAGSTAT {
    container 'ghcr.io/bf528/samtools:latest'
    publishDir params.outdir, mode: "copy", pattern: '*_flagstat.txt'
	label 'process_medium'

    input:
    tuple val(sample), val(name), path(bam)
    output:
    tuple val(sample), val(name), path("${sample}_flagstat.txt"), emit: flagstat

    script:
    """
    samtools flagstat ${bam} > ${sample}_flagstat.txt
    """
    stub:
    """
    touch ${sample}_flagstat.txt
    """
}
