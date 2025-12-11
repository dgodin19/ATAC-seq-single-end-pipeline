#!/usr/bin/env nextflow

process PREPARE_DIFFBIND_SAMPLESHEET {
    label 'process_single'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment),
        path(filtered_bam), path(filtered_bam_index), path(narrowPeak),
        val(factor), val(cell_type), val(ko_wt), val(replicate)

    output:
    path "*.csv"

    script:
    // Unique filename to avoid collisions
    def file_id = "${cell_type}_${ko_wt}_${replicate}"

    // Condition WITHOUT replicate (correct for DiffBind)
    def diffbind_condition = "${cell_type}_${ko_wt}"

    """
    echo "Generating sample sheet for ${file_id}"
    echo "SampleID,Tissue,Factor,Condition,Replicate,bamReads,Peaks,PeakCaller" \
        > ${file_id}_diffbind_samplesheet.csv

    echo "${run},cells,${factor},${diffbind_condition},${replicate}," \
         "${filtered_bam.toRealPath()},${narrowPeak.toRealPath()},narrowPeak" \
         >> ${file_id}_diffbind_samplesheet.csv

    cat ${file_id}_diffbind_samplesheet.csv
    """
}