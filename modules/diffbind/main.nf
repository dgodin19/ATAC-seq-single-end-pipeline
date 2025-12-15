#!/usr/bin/env nextflow

process DIFFBIND {
    label 'process_single'
    conda 'envs/diffbind_env.yml'
    publishDir params.outdir, mode: 'copy', pattern: '*'

    input:
    path sample_sheet

    output:
    path "*_significant_peaks_p001.bed", emit: significant_peaks
    tuple path("*_gain.bed"), path("*_loss.bed"), emit: gain_loss
    path "*_diffbind_results.csv", emit: results_csv

    script:
    def output_prefix = sample_sheet.baseName
    """
    echo "Processing sample sheet: ${sample_sheet}"
    ls -l ${sample_sheet}
    cat ${sample_sheet}
    diffbind.R ${sample_sheet} ${output_prefix} DESEQ2
    """
}
