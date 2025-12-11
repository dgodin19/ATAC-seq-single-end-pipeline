#!/usr/bin/env nextflow

process MERGE_SAMPLESHEETS {
  publishDir params.outdir, mode: 'copy'

  input:
    tuple val(cell_type), path(csvs)

  output:
    path "${cell_type}_diffbind_samplesheet.csv"

  script:
    """
    head -n1 \$(ls $csvs | head -n1) > ${cell_type}_diffbind_samplesheet.csv
    for f in $csvs; do tail -n +2 "\$f"; done >> ${cell_type}_diffbind_samplesheet.csv
    """
}
