#!/usr/bin/env nextflow

process FRIP_CALC {
    label 'process_medium'
    conda 'envs/bedtools_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '*.txt'

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(filtered_bam), 
    path(filtered_bam_index),
    path(peaks_bed),
    path(summits_bed),
    path(peaks_xls)

    output:
    path("${run}_frip.txt")

    script:
    """
    total_reads=\$(samtools view -c ${filtered_bam})
    in_peaks=\$(bedtools intersect -u -a ${filtered_bam} -b ${peaks_bed} | wc -l)
    frip=\$(echo "scale=4; \$in_peaks / \$total_reads" | bc)

    echo -e "run,frip" > ${run}_frip.txt
    echo -e "${run},\$frip" >> ${run}_frip.txt
    """
}
