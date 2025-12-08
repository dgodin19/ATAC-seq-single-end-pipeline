#!/usr/bin/env nextflow

process REMOVE_BLACKLIST {
    label 'process_medium'
    conda 'envs/bedtools_env.yml'

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment),
          path(narrowPeak),
          path(summits_bed),
          path(xls)
    path blacklist

    output:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path("${run}.blacklist_removed.bed")

    script:
    """
    bedtools intersect -v -a ${summits_bed} -b ${blacklist} > ${run}.blacklist_removed.bed

    # Sort & index the new BAM
    # samtools sort -o ${run}.blacklist_removed.sorted.bam ${run}.blacklist_removed.bam
    # samtools index ${run}.blacklist_removed.sorted.bam
    # rm ${run}.blacklist_removed.bam
    # samtools view -c ${run}.blacklist_removed.sorted.bam
    """

    stub:
    """
    touch ${run}.blacklist_removed.bed
    """
}
