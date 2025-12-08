#!/usr/bin/env nextflow


process MERGE_PEAKS {
    label 'process_single'
    conda 'envs/bedtools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.bed'

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(shifted_blacklist_removed_bam),  path(shifted_blacklist_removed_bam_index)

    output:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment),
          path("${run}_peaks.narrowPeak"),
          path("${run}_summits.bed"),
          path("${run}_peaks.xls"),
          path("${run}_treat_pileup.bdg")

    script:
    """
    samtools view -c ${shifted_blacklist_removed_bam}
    macs2 callpeak \
        -t ${shifted_blacklist_removed_bam} \
        -f BAM \
        -g mm \
        --nomodel \
        --keep-dup auto \
        --extsize 147 \
        -n ${run}
    """
}
