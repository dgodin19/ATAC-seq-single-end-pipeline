#!/usr/bin/env nextflow

process ANNOTATE {
    conda 'envs/homer_env.yml'
    publishDir "${params.outdir}/annotation/${condition}", mode: "copy"
    label 'process_high'

    input:
    path filtered_bed
    path genome
    path gtf
    val condition

    output:
    path("annotated_${condition}_peaks.txt")
    path("annotate_${condition}_log.txt")

    script:
    """
    # Use condition in output filenames
    awk 'BEGIN{OFS="\t"} {print "Peak_"NR, \$1, \$2, \$3, \$5, \$6}' ${filtered_bed} > converted_${condition}_peaks.txt

    annotatePeaks.pl ${filtered_bed} ${genome} \
    -gtf ${gtf} > annotated_${condition}_peaks.txt 2> annotate_${condition}_log.txt

    # Prepend the correct header
    (echo -e "PeakID\tChr\tStart\tEnd\tStrand\tPeak Score\tAnnotation\tGene Name\tGene Alias\tDistance to TSS\tNearest PromoterID\tEntrez ID\tNearest Unigene\tNearest Refseq\tNearest Ensembl\tGene Description\tGene Type"; 
     cat annotated_${condition}_peaks.txt) > temp_peaks.txt

    # Replace the original file
    mv temp_peaks.txt annotated_${condition}_peaks.txt

    # Check output
    echo "Annotation output:"
    head -n 5 annotated_${condition}_peaks.txt
    wc -l annotated_${condition}_peaks.txt

    # Show any errors or logs
    cat annotate_${condition}_log.txt
    """

    stub:
    """
    echo -e "PeakID\tChr\tStart\tEnd\tStrand\tAnnotation\tGene Name\tGene Alias\tGene Description" > annotated_${condition}_peaks.txt
    echo -e "Peak_1\tchr1\t1000\t2000\t+\tPromoter\tGene1\tAlias1\tDescription1" >> annotated_${condition}_peaks.txt
    touch annotate_${condition}_log.txt
    """
}