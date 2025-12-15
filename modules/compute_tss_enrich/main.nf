#!/usr/bin/env nextflow

process COMPUTE_TSS_ENRICH {
    conda 'envs/pandas_env.yml'
    label 'process_low'
    publishDir params.outdir, mode: "copy", pattern: '*.txt'

    input:
    tuple val (celltype), path (signal_coverage_tab)

    output:
    path("${celltype}_tss_enrichment.txt")

    script:
    """
    compute_tss_enrichment.py \
       --input ${signal_coverage_tab} \
       --output ${celltype}_tss_enrichment.txt
    """
}
