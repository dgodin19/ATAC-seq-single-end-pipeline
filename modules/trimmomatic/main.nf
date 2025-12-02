#!/usr/bin/env nextflow

process TRIM {
    label 'process_low'
    conda 'envs/trimmomatic_env.yml'
	publishDir params.outdir, mode: "copy", pattern: '*.log'
	
    input:
	tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(fastq)
    path(adapters)

	output:
	tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path('*_trim.log'), emit: log
	tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), emit: trim

    
    script:
    """
    trimmomatic SE -phred33 ${fastq} ${run}_trimmed.fastq.gz \
        ILLUMINACLIP:${adapters}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        2> ${run}_trim.log
    """ 
    
    stub:
    """
    touch ${run}_stub_trim.log
    touch ${run}_stub_trimmed.fastq.gz
    """
}
