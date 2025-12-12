#!/usr/bin/env nextflow

process PLOTPROFILE {
    conda 'envs/deeptools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.png'
	label 'process_medium'
    
    input:
    path(compute_matrix)

    output:
    path('signal_coverage.png'), emit: signal_png
    path('signal_coverage.tab'), emit: signal_tab

    script:
    """
    plotProfile \
    -m $compute_matrix \
    -o signal_coverage.png \
    --outFileNameData signal_coverage.tab
    """


    stub:
    """
    touch signal_coverage.png
    """
}