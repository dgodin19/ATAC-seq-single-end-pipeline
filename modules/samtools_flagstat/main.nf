#!/usr/bin/env nextflow

process SAMTOOLS_FLAGSTAT {
    label 'process_medium'
    conda 'envs/samtools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*_flagstat.txt'
	
    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(bam)
    output:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path("${run}_flagstat.txt"), emit: flagstat

    script:
    """
    samtools flagstat ${bam} > ${run}_flagstat.txt
    """
    stub:
    """
    touch ${run}_flagstat.txt
    """
}
