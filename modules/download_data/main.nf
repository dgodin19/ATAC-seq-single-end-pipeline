#!/usr/bin/env nextflow

process DOWNLOAD_SAMPLES {
    label 'process_single'
    conda 'envs/ncbidatasets_env.yml'
    input:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment)

    output:
    tuple val(run), val(biosample), val(samplename), val(library_layout), val(library_source), val(experiment), path("${run}*_R*.fastq")

    script:
    """
    # Download as before, using run, library_layout, etc.
    prefetch ${run}
    if [ "${library_layout}" == "PAIRED" ]; then
        fasterq-dump --split-files --outdir . ${run}
        mv ${run}_1.fastq ${run}_R1.fastq
        mv ${run}_2.fastq ${run}_R2.fastq
    else
        fasterq-dump --outdir . ${run}
        mv ${run}.fastq ${run}_R1.fastq
    fi
    rm -rf ${run}
    """
}