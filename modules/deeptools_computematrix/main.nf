#!/usr/bin/env nextflow

process COMPUTEMATRIX {
    conda 'envs/deeptools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.gz'
    label 'process_medium'

    input:
    val(celltype) 
    path(bigwigs) 
    path(sigpeak_bed)

    output:
    tuple val(celltype), path("${celltype}_matrix.gz")

    script:
    """
    computeMatrix reference-point \
      --referencePoint center \
      -S ${bigwigs.join(' ')} \
      -R $sigpeak_bed \
      -b 500 -a 500 \
      -o ${celltype}_matrix.gz
    """

    stub:
    """
    touch ${celltype}_matrix.gz
    """
}