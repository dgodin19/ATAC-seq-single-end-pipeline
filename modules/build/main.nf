#!/usr/bin/env nextflow

process BOWTIE2_BUILD {
    label 'process_high'
    conda 'envs/bowtie2_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.html'
	
    input:
    tuple val(name), path(genome)

    output:
    tuple val(name), path('bowtie2_index/'), emit: index

    script:
    """
    mkdir -p bowtie2_index
    bowtie2-build ${genome} bowtie2_index/${name}
    """
    stub:
    """
    mkdir bowtie2_index
    touch bowtie2_index/stub.bt2*
    """
}
