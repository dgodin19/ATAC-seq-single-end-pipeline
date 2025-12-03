#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {
    label 'process_high'
    conda 'envs/bowtie2_env.yml'

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path(fastq)
    tuple val(name), path(index)

    output:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path("${run}.filtered.bam"), path("${run}.filtered.bam.bai"), emit: filtered_bam
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment),  path("${run}.bam"), emit: samtools_flagstat

    script:
    """
    bowtie2 \
    -x ${index}/${name} \
    -U ${fastq} \
    -S ${run}.sam \
    --very-sensitive \
    -p ${task.cpus} \
    -q
    
    # Remove -m 1 and --maxins 2000 flags
    # Add threads with -p
    # Consider adding --no-unal to reduce output size

    # Convert to BAM and filter
    samtools view -q 30 -F 4 -b ${run}.sam > ${run}.bam
    samtools sort ${run}.bam -o ${run}.sorted.bam
    rm ${run}.sam

    samtools view -h ${run}.sorted.bam | \
        grep -v "chrM" | \
        samtools view -bS > ${run}.filtered.bam
    samtools index ${run}.filtered.bam
    rm ${run}.sorted.bam
    """
    stub:
    """
    touch ${run}.bam
    touch ${run}.filtered.bam
    """
}