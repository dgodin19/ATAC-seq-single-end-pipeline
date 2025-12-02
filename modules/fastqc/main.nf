#!/usr/bin/env nextflow

process FASTQC {
    label 'process_high'
    conda 'envs/fastqc_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.{html,zip}'
    
    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(fastq)

    output:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path('*_fastqc.html'), emit: html
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path('*_fastqc.zip'), emit: zip

    script:
    """
    fastqc $fastq -t $task.cpus
    """

    stub:
    """
    touch ${run}_fastqc.html
    touch ${run}_fastqc.zip
    """
}


