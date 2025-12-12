#!/usr/bin/env nextflow

process BAMCOVERAGE {
    conda 'envs/deeptools_env.yml'
	label 'process_medium'

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(filtered_bam), path(filtered_bam_index)

    output:
    tuple val(run), path("${run}.bw")

    script:
    """
    bamCoverage \
        -b ${filtered_bam} \
        -o ${run}.bw \
        --binSize 10 \
        --normalizeUsing CPM
    """

    stub:
    """
    touch ${run}.bw
    """
}