#!/usr/bin/env nextflow

process SHIFT_ALIGN {
    label 'process_high'
    conda 'envs/deeptools_env.yml'

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(filtered_bam), path(filtered_bam_index)

    output:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path("${run}_shifted.bam"), path("${run}_shifted.bam.bai")

    script:
    """
    alignmentSieve -b ${filtered_bam} -o ${run}_shifted.bam  --ATACshift --minMappingQuality 10 --filterMetrics log.txt
    samtools index ${run}_shifted.bam
    """
}