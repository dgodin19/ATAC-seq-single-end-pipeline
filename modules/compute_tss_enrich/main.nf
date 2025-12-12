#!/usr/bin/env nextflow

process COMPUTE_TSS_ENRICH {
    conda 'envs/pandas_env.yml'
    label 'process_low'
    publishDir params.outdir, mode: "copy", pattern: '*'

    input:
    path signal_coverage

    output:
    path 'tss_enrichment.txt'

    script:
    """
    ./compute_tss_enrichment.py --input ${signal_coverage} --output tss_enrichment.txt
    """
}
