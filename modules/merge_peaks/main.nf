#!/usr/bin/env nextflow


process MERGE_PEAKS {
    label 'process_merge'
    conda 'envs/bedtools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.bed'

    input:
    path bed_files
    val condition

    output:
    path "merged_${condition}_peaks.bed"

    script:
    """
    cat ${bed_files} \
        | sort -k1,1 -k2,2n \
        | bedtools merge > merged_${condition}_peaks.bed
    """
}