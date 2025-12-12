#!/usr/bin/env nextflow

process ANNOTATE {
    input:
    path filtered_bed
    path genome
    path gtf

    output:
    path("annotated_*_peaks.txt")

    script:
    // derive condition from filename
    def condition = filtered_bed.baseName.replace('_significant_peaks_p001','')

    """
    annotatePeaks.pl ${filtered_bed} ${genome} -gtf ${gtf} > annotated_${condition}_peaks.txt
    """
}