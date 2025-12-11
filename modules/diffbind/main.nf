#!/usr/bin/env nextflow

process DIFFBIND {
    label 'process_single'
    conda 'envs/diffbind_env.yml'
    publishDir params.outdir, mode: 'copy'

    input:
    path sample_sheet

    output:
    path "*_significant_peaks_p001.bed"

    script:
    def output_prefix = sample_sheet.baseName
    """
    echo "Processing sample sheet: ${sample_sheet}"
    ls -l ${sample_sheet}
    cat ${sample_sheet}
    diffbind.R ${sample_sheet} ${output_prefix} DESEQ2
    """
}
