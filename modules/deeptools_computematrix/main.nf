#!/usr/bin/env nextflow

process COMPUTEMATRIX {
    conda 'envs/deeptools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.gz'
	label 'process_medium'

    input:
    val(bigwigs) 
    path(bed) 
    val(window)

    output:
    path('matrix.gz')

    script:

    """
    computeMatrix reference-point \
    -S ${bigwigs.join(' ')} \
    -R ${bed} \
    --referencePoint TSS \
    -b ${window} \
    -a ${window} \
    -o matrix.gz
    """

    stub:
    """
    touch matrix.gz
    """
}