#!/usr/bin/env nextflow

process SAMTOOLS_GENOME {
    conda 'envs/samtools_env.yml'
    label 'process_medium'
    
    input:
    path genome

    output:
    tuple path(genome), path("${genome}.fai")

    script:
    """
    samtools faidx ${genome}
    """
}
