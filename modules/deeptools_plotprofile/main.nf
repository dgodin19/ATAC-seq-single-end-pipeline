#!/usr/bin/env nextflow

process PLOTPROFILE {
    conda 'envs/deeptools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*'
	label 'process_medium'
    
    input:
    tuple val(celltype), path(compute_matrix)

    output:
    path("${celltype}_signal_coverage.png"), emit: signal_png
    tuple val(celltype), path("${celltype}_signal_coverage.tab"), emit: signal_tab

    script:
    """
    plotProfile \
      -m ${compute_matrix} \
      -o ${celltype}_signal_coverage.png \
      --outFileNameData ${celltype}_signal_coverage.tab
    """


    stub:
    """
    touch signal_coverage.png
    """
}