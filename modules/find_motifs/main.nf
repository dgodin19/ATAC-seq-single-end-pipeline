#!/usr/bin/env nextflow

process FIND_MOTIFS {
    conda 'envs/homer_env.yml'
    publishDir params.outdir, mode: "copy", pattern: '{motif_output_*,findMotifs_*.log}'
    label 'process_high'

    input:
    path filtered_bed
    tuple path(genome), path(genome_idx)

    output:
    path("motif_output_*")
    path("findMotifs_*.log")

    script:
    // derive condition from filename
    def condition = filtered_bed.baseName.replace('_significant_peaks_p001','')

    """
    # Create output directory for this condition
    mkdir -p motif_output_${condition}
    awk '{OFS="\t"; print \$1,\$2,\$3,".",\$5,\$6}' ${filtered_bed} > ${condition}.bed

    # Run HOMER
    findMotifsGenome.pl ${condition}.bed ${genome} motif_output_${condition} \
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