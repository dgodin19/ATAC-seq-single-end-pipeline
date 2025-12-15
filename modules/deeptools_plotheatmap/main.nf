#!/usr/bin/env nextflow

process PLOTHEATMAP {
    label 'process_medium'
    conda 'envs/deeptools_env.yml'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(celltype), path(matrix)

    output:
    path("*.png")

    script:
    """
    plotHeatmap -m $matrix --kmeans 3 --colorMap Blues -o ${celltype}_diff_peaks_heatmap.png 
    """

    stub:
    """
    touch stub_diff_peaks_heatmap.png
    """
}
