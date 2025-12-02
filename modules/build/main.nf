#!/usr/bin/env nextflow

process BOWTIE2_BUILD {
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir params.outdir, mode: "copy", pattern: '*.html'
	label 'process_high'

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
