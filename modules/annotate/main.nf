#!/usr/bin/env nextflow

process ANNOTATE {
    conda 'envs/homer_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.txt'
    
    input:
    path filtered_bed
    tuple path(genome), path(genome_idx)
    path gtf

    output:
    path("annotated_*_peaks.txt")

    script:
    // derive condition from filename
    def condition = filtered_bed.baseName.replace('_significant_peaks_p001','')

    """
    awk '{OFS="\t"; print \$1,\$2,\$3,".",\$5,\$6}' ${filtered_bed} > ${condition}.bed
    annotatePeaks.pl ${condition}.bed ${genome} -gtf ${gtf} > annotated_${condition}_peaks.txt
    """
}