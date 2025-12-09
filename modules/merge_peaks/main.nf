#!/usr/bin/env nextflow


process MERGE_PEAKS {
    label 'process_single'
    conda 'envs/bedtools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.bed'

    input:
    path bed_list

    output:
    path "consensus_peaks.bed"

    script:
    """
    cat *.bed \
    | sort -k1,1 -k2,2n \
    | bedtools merge \
    > consensus_peaks.bed
    """
}
