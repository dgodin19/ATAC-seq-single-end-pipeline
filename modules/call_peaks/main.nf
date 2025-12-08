#!/usr/bin/env nextflow


process CALL_PEAKS {
    label 'process_single'
    conda 'envs/macs2_env.yml'
    publishDir params.outdir, mode: "copy", pattern: ['*.narrowPeak','*.bed','*.xls']

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(filtered_bam),  path(filtered_bam_index)

    output:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment),
          path("${run}_peaks.narrowPeak"),
          path("${run}_summits.bed"),
          path("${run}_peaks.xls")

    script:
    """
    samtools view -c ${filtered_bam}
    macs2 callpeak \
        -t ${filtered_bam} \
        -f BAM \
        -g mm \
        --nomodel \
        --keep-dup auto \
        --extsize 147 \
        -n ${run}
    """
}