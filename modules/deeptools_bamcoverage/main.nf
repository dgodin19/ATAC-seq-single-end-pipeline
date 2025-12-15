#!/usr/bin/env nextflow

process BAMCOVERAGE {
    conda 'envs/deeptools_env.yml'
	label 'process_medium'

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(filtered_bam), path(filtered_bam_index)
    val run2cond_map 

    output:
    tuple val(run), path("${run2cond_map[run]}.bw")

    script:
    """
    bamCoverage \
        -b ${filtered_bam} \
        -o ${run2cond_map[run]}.bw \
        --binSize 10 \
        --normalizeUsing CPM
    """

    stub:
    """
    touch ${run2cond[run]}.bw
    """
}