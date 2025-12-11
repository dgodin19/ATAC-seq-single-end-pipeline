#!/usr/bin/env nextflow

process FIND_MOTIFS {
    conda 'envs/homer_env.yml'
    publishDir "${params.outdir}/motifs/${condition}", mode: "copy"
    label 'process_high'
    
    input:
    path filtered_bed
    path genome
    val condition

    output:
    path("motif_output_${condition}")
    path("findMotifs_${condition}.log")

    script:
    """
    mkdir -p motif_output_${condition}

    findMotifsGenome.pl ${filtered_bed} ${genome} motif_output_${condition} \
        -size 200 \
        -len 8,10,12 \
        -preparse \
        -p ${task.cpus} 2>&1 | tee findMotifs_${condition}.log
    """

    stub:
    """
    mkdir -p motif_output_${condition}/homerResults
    mkdir -p motif_output_${condition}/knownResults
    touch motif_output_${condition}/homerMotifs.all.motifs
    touch motif_output_${condition}/homerResults.html
    touch findMotifs_${condition}.log
    """
}


