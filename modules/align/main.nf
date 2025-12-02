#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir params.outdir, mode: "copy", pattern: '*.bam'
    label 'process_medium'

    input:
    tuple val(sample), path(fastq)
    tuple val(name), path(index)

    output:
    tuple val(sample), val(name), path("${sample}.bam"), emit: bam

    script:
    """
    bowtie2 \
    -x ${index}/${name} \
    -U ${fastq} \
    -S ${sample}.sam \
    --very-sensitive \
    --no-discordant \
    --no-mixed \
    -m 1 \
    --maxins 2000 \
    -q
    
    # Convert to BAM and filter
    samtools view -q 30 -F 4 -f 2 -b ${sample}.sam > ${sample}.bam
    samtools sort ${sample}.bam -o ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam
    rm ${sample}.sam

    # Optional: additional ATAC-seq specific filtering
    samtools view -h ${sample}.sorted.bam | \
    grep -v "chrM" | \
    samtools view -bS > ${sample}.filtered.bam
    """
   
    stub:
    """
    touch ${sample}.bam
    """
}